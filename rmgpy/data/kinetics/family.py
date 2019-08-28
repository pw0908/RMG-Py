#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This module contains functionality for working with kinetics families.
"""
from __future__ import division

import codecs
import itertools
import logging
import multiprocessing as mp
import os.path
import random
import re
import warnings
from collections import OrderedDict
from copy import deepcopy

import numpy as np
from sklearn.model_selection import KFold

from rmgpy import settings
from rmgpy.constraints import failsSpeciesConstraints
from rmgpy.data.base import Database, Entry, LogicNode, LogicOr, ForbiddenStructures, getAllCombinations
from rmgpy.data.kinetics.common import saveEntry, find_degenerate_reactions, generate_molecule_combos, \
                                       ensure_independent_atom_ids
from rmgpy.data.kinetics.depository import KineticsDepository
from rmgpy.data.kinetics.groups import KineticsGroups
from rmgpy.data.kinetics.rules import KineticsRules
from rmgpy.exceptions import ActionError, DatabaseError, InvalidActionError, KekulizationError, KineticsError, \
                             ForbiddenStructureException, UndeterminableKineticsError
from rmgpy.kinetics import Arrhenius, SurfaceArrhenius, SurfaceArrheniusBEP, StickingCoefficient, \
                           StickingCoefficientBEP, ArrheniusBM
from rmgpy.kinetics.uncertainties import RateUncertainty, rank_accuracy_map
from rmgpy.molecule import Bond, GroupBond, Group, Molecule
from rmgpy.molecule.atomtype import atomTypes
from rmgpy.reaction import Reaction, same_species_lists
from rmgpy.species import Species


################################################################################

class TemplateReaction(Reaction):
    """
    A Reaction object generated from a reaction family template. In addition
    to attributes inherited from :class:`Reaction`, this class includes the
    following attributes:

    ============ ========================= =====================================
    Attribute    Type                      Description
    ============ ========================= =====================================
    `family`     ``str``                   The kinetics family that the reaction was created from.
    `estimator`  ``str``                   Whether the kinetics came from rate rules or group additivity.
    `reverse`    :class:`TemplateReaction` The reverse reaction, for families that are their own reverse.
    `is_forward`  ``bool``                 Whether the reaction was generated in the forward direction of the family.
    ============ ========================= =====================================
    """

    def __init__(self,
                 index=-1,
                 reactants=None,
                 products=None,
                 specificCollider=None,
                 kinetics=None,
                 reversible=True,
                 transitionState=None,
                 duplicate=False,
                 degeneracy=1,
                 pairs=None,
                 family=None,
                 template=None,
                 estimator=None,
                 reverse=None,
                 is_forward=None,
                 ):
        Reaction.__init__(self,
                          index=index,
                          reactants=reactants,
                          products=products,
                          specificCollider=specificCollider,
                          kinetics=kinetics,
                          reversible=reversible,
                          transitionState=transitionState,
                          duplicate=duplicate,
                          degeneracy=degeneracy,
                          pairs=pairs,
                          is_forward=is_forward,
                          )
        self.family = family
        self.template = template
        self.estimator = estimator
        self.reverse = reverse

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (TemplateReaction, (self.index,
                                   self.reactants,
                                   self.products,
                                   self.specificCollider,
                                   self.kinetics,
                                   self.reversible,
                                   self.transitionState,
                                   self.duplicate,
                                   self.degeneracy,
                                   self.pairs,
                                   self.family,
                                   self.template,
                                   self.estimator,
                                   self.reverse,
                                   self.is_forward
                                   ))

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'TemplateReaction('
        if self.index != -1: string += 'index={0:d}, '.format(self.index)
        if self.label != '': string += 'label={0!r}, '.format(self.label)
        if self.reactants is not None: string += 'reactants={0!r}, '.format(self.reactants)
        if self.products is not None: string += 'products={0!r}, '.format(self.products)
        if self.specificCollider is not None: string += 'specificCollider={0!r}, '.format(self.specificCollider)
        if self.kinetics is not None: string += 'kinetics={0!r}, '.format(self.kinetics)
        if not self.reversible: string += 'reversible={0}, '.format(self.reversible)
        if self.transitionState is not None: string += 'transitionState={0!r}, '.format(self.transitionState)
        if self.duplicate: string += 'duplicate={0}, '.format(self.duplicate)
        if self.degeneracy != 1: string += 'degeneracy={0:.1f}, '.format(self.degeneracy)
        if self.pairs is not None: string += 'pairs={0}, '.format(self.pairs)
        if self.family: string += "family='{}', ".format(self.family)
        if self.template: string += "template={}, ".format(self.template)
        if self.comment != '': string += 'comment={0!r}, '.format(self.comment)
        string = string[:-2] + ')'
        return string

    def getSource(self):
        """
        Return the database that was the source of this reaction. For a
        TemplateReaction this should be a KineticsGroups object.
        """
        return self.family

    def copy(self):
        """
        creates a new instance of TemplateReaction
        """
        other = TemplateReaction.__new__(TemplateReaction)

        # this was copied from Reaction.copy class
        other.index = self.index
        other.label = self.label
        other.reactants = []
        for reactant in self.reactants:
            other.reactants.append(reactant.copy(deep=True))
        other.products = []
        for product in self.products:
            other.products.append(product.copy(deep=True))
        other.specificCollider = self.specificCollider
        other.degeneracy = self.degeneracy
        other.kinetics = deepcopy(self.kinetics)
        other.reversible = self.reversible
        other.transitionState = deepcopy(self.transitionState)
        other.duplicate = self.duplicate
        other.pairs = deepcopy(self.pairs)

        # added for TemplateReaction information
        other.family = self.family
        other.template = self.template
        other.estimator = self.estimator
        other.reverse = self.reverse
        other.is_forward = self.is_forward

        return other


################################################################################

class ReactionRecipe(object):
    """
    Represent a list of actions that, when executed, result in the conversion
    of a set of reactants to a set of products. There are currently five such
    actions:

    ============= ============================= ================================
    Action Name   Arguments                     Description
    ============= ============================= ================================
    CHANGE_BOND   `center1`, `order`, `center2` change the bond order of the bond between `center1` and `center2` by `order`; do not break or form bonds
    FORM_BOND     `center1`, `order`, `center2` form a new bond between `center1` and `center2` of type `order`
    BREAK_BOND    `center1`, `order`, `center2` break the bond between `center1` and `center2`, which should be of type `order`
    GAIN_RADICAL  `center`, `radical`           increase the number of free electrons on `center` by `radical`
    LOSE_RADICAL  `center`, `radical`           decrease the number of free electrons on `center` by `radical`
    GAIN_PAIR     `center`, `pair`              increase the number of lone electron pairs on `center` by `pair`
    LOSE_PAIR     `center`, `pair`              decrease the number of lone electron pairs on `center` by `pair`
    ============= ============================= ================================

    The actions are stored as a list in the `actions` attribute. Each action is
    a list of items; the first is the action name, while the rest are the
    action parameters as indicated above.
    """

    def __init__(self, actions=None):
        self.actions = actions or []

    def addAction(self, action):
        """
        Add an `action` to the reaction recipe, where `action` is a list
        containing the action name and the required parameters, as indicated in
        the table above.
        """
        self.actions.append(action)

    def getReverse(self):
        """
        Generate a reaction recipe that, when applied, does the opposite of
        what the current recipe does, i.e., it is the recipe for the reverse
        of the reaction that this is the recipe for.
        """
        other = ReactionRecipe()
        for action in self.actions:
            if action[0] == 'CHANGE_BOND':
                other.addAction(['CHANGE_BOND', action[1], str(-int(action[2])), action[3]])
            elif action[0] == 'FORM_BOND':
                other.addAction(['BREAK_BOND', action[1], action[2], action[3]])
            elif action[0] == 'BREAK_BOND':
                other.addAction(['FORM_BOND', action[1], action[2], action[3]])
            elif action[0] == 'LOSE_RADICAL':
                other.addAction(['GAIN_RADICAL', action[1], action[2]])
            elif action[0] == 'GAIN_RADICAL':
                other.addAction(['LOSE_RADICAL', action[1], action[2]])
            elif action[0] == 'LOSE_PAIR':
                other.addAction(['GAIN_PAIR', action[1], action[2]])
            elif action[0] == 'GAIN_PAIR':
                other.addAction(['LOSE_PAIR', action[1], action[2]])
        return other

    def __apply(self, struct, doForward, unique):
        """
        Apply the reaction recipe to the set of molecules contained in
        `structure`, a single Structure object that contains one or more
        structures. The `doForward` parameter is used to indicate
        whether the forward or reverse recipe should be applied. The atoms in
        the structure should be labeled with the appropriate atom centers.
        """

        pattern = isinstance(struct, Group)
        struct.props['validAromatic'] = True

        for action in self.actions:
            if action[0] in ['CHANGE_BOND', 'FORM_BOND', 'BREAK_BOND']:

                # We are about to change the connectivity of the atoms in
                # struct, which invalidates any existing vertex connectivity
                # information; thus we reset it
                struct.resetConnectivityValues()

                label1, info, label2 = action[1:]

                if label1 != label2:
                    # Find associated atoms
                    atom1 = struct.getLabeledAtom(label1)[0]
                    atom2 = struct.getLabeledAtom(label2)[0]
                else:
                    atoms = struct.getLabeledAtom(label1)  # should never have more than two if this action is valid
                    if len(atoms) > 2:
                        raise InvalidActionError('Invalid atom labels encountered.')
                    atom1, atom2 = atoms

                if atom1 is None or atom2 is None or atom1 is atom2:
                    raise InvalidActionError('Invalid atom labels encountered.')

                # Apply the action
                if action[0] == 'CHANGE_BOND':
                    info = int(info)
                    bond = struct.getBond(atom1, atom2)
                    if bond.isBenzene():
                        struct.props['validAromatic'] = False
                    if doForward:
                        atom1.applyAction(['CHANGE_BOND', label1, info, label2])
                        atom2.applyAction(['CHANGE_BOND', label1, info, label2])
                        bond.applyAction(['CHANGE_BOND', label1, info, label2])
                    else:
                        atom1.applyAction(['CHANGE_BOND', label1, -info, label2])
                        atom2.applyAction(['CHANGE_BOND', label1, -info, label2])
                        bond.applyAction(['CHANGE_BOND', label1, -info, label2])
                elif (action[0] == 'FORM_BOND' and doForward) or (action[0] == 'BREAK_BOND' and not doForward):
                    if struct.hasBond(atom1, atom2):
                        raise InvalidActionError('Attempted to create an existing bond.')
                    if info not in (1, 0):  # Can only form single or vdW bonds
                        raise InvalidActionError('Attempted to create bond of type {:!r}'.format(info))
                    bond = GroupBond(atom1, atom2, order=[info]) if pattern else Bond(atom1, atom2, order=info)
                    struct.addBond(bond)
                    atom1.applyAction(['FORM_BOND', label1, info, label2])
                    atom2.applyAction(['FORM_BOND', label1, info, label2])
                elif (action[0] == 'BREAK_BOND' and doForward) or (action[0] == 'FORM_BOND' and not doForward):
                    if not struct.hasBond(atom1, atom2):
                        raise InvalidActionError('Attempted to remove a nonexistent bond.')
                    bond = struct.getBond(atom1, atom2)
                    struct.removeBond(bond)
                    atom1.applyAction(['BREAK_BOND', label1, info, label2])
                    atom2.applyAction(['BREAK_BOND', label1, info, label2])

            elif action[0] in ['LOSE_RADICAL', 'GAIN_RADICAL']:

                label, change = action[1:]
                change = int(change)

                # Find associated atom
                atoms = struct.getLabeledAtom(label)
                for atom in atoms:
                    if atom is None:
                        raise InvalidActionError('Unable to find atom with label "{0}" while applying '
                                                 'reaction recipe.'.format(label))

                    # Apply the action
                    for i in range(change):
                        if (action[0] == 'GAIN_RADICAL' and doForward) or (action[0] == 'LOSE_RADICAL' and not doForward):
                            atom.applyAction(['GAIN_RADICAL', label, 1])
                        elif (action[0] == 'LOSE_RADICAL' and doForward) or (action[0] == 'GAIN_RADICAL' and not doForward):
                            atom.applyAction(['LOSE_RADICAL', label, 1])

            elif action[0] in ['LOSE_PAIR', 'GAIN_PAIR']:

                label, change = action[1:]
                change = int(change)

                # Find associated atom
                atoms = struct.getLabeledAtom(label)

                for atom in atoms:
                    if atom is None:
                        raise InvalidActionError('Unable to find atom with label "{0}" while applying '
                                                 'reaction recipe.'.format(label))

                    # Apply the action
                    for i in range(change):
                        if (action[0] == 'GAIN_PAIR' and doForward) or (action[0] == 'LOSE_PAIR' and not doForward):
                            atom.applyAction(['GAIN_PAIR', label, 1])
                        elif (action[0] == 'LOSE_PAIR' and doForward) or (action[0] == 'GAIN_PAIR' and not doForward):
                            atom.applyAction(['LOSE_PAIR', label, 1])

            else:
                raise InvalidActionError('Unknown action "' + action[0] + '" encountered.')

    def applyForward(self, struct, unique=True):
        """
        Apply the forward reaction recipe to `molecule`, a single
        :class:`Molecule` object.
        """
        return self.__apply(struct, True, unique)

    def applyReverse(self, struct, unique=True):
        """
        Apply the reverse reaction recipe to `molecule`, a single
        :class:`Molecule` object.
        """
        return self.__apply(struct, False, unique)


################################################################################


class KineticsFamily(Database):
    """
    A class for working with an RMG kinetics family: a set of reactions with 
    similar chemistry, and therefore similar reaction rates. The attributes 
    are:

    =================== =============================== ========================
    Attribute           Type                            Description
    =================== =============================== ========================
    `reverse`           ``string``                      The name of the reverse reaction family
    `reversible`        `Boolean`                       Is family reversible? (True by default)
    `forwardTemplate`   :class:`Reaction`               The forward reaction template
    `forwardRecipe`     :class:`ReactionRecipe`         The steps to take when applying the forward reaction to a set of reactants
    `reverseTemplate`   :class:`Reaction`               The reverse reaction template
    `reverseRecipe`     :class:`ReactionRecipe`         The steps to take when applying the reverse reaction to a set of reactants
    `forbidden`         :class:`ForbiddenStructures`    (Optional) Forbidden product structures in either direction
    `ownReverse`        `Boolean`                       It's its own reverse?
    'boundaryAtoms'     list                            Labels which define the boundaries of end groups in backbone/end families
    `treeDistances`     dict                            The default distance from parent along each tree, if not set default is 1 for every tree
    ------------------- ------------------------------- ------------------------
    `groups`            :class:`KineticsGroups`         The set of kinetics group additivity values
    `rules`             :class:`KineticsRules`          The set of kinetics rate rules from RMG-Java
    `depositories`      ``list``                        A set of additional depositories used to store kinetics data from various sources
    =================== =============================== ========================

    There are a few reaction families that are their own reverse (hydrogen
    abstraction and intramolecular hydrogen migration); for these
    `reverseTemplate` and `reverseRecipe` will both be ``None``.
    """

    def __init__(self,
                 entries=None,
                 top=None,
                 label='',
                 name='',
                 reverse='',
                 reversible=True,
                 shortDesc='',
                 longDesc='',
                 forwardTemplate=None,
                 forwardRecipe=None,
                 reverseTemplate=None,
                 reverseRecipe=None,
                 forbidden=None,
                 boundaryAtoms=None,
                 treeDistances=None
                 ):
        Database.__init__(self, entries, top, label, name, shortDesc, longDesc)
        self.reverse = reverse
        self.reversible = reversible
        self.forwardTemplate = forwardTemplate
        self.forwardRecipe = forwardRecipe
        self.reverseTemplate = reverseTemplate
        self.reverseRecipe = reverseRecipe
        self.forbidden = forbidden
        self.ownReverse = forwardTemplate is not None and reverseTemplate is None
        self.boundaryAtoms = boundaryAtoms
        self.treeDistances = treeDistances

        # Kinetics depositories of training and test data
        self.groups = None
        self.rules = None
        self.depositories = []

    def __repr__(self):
        return '<ReactionFamily "{0}">'.format(self.label)

    def loadOld(self, path):
        """
        Load an old-style RMG kinetics group additivity database from the
        location `path`.
        """
        warnings.warn("The old kinetics databases are no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        self.label = os.path.basename(path)
        self.name = self.label

        self.groups = KineticsGroups(label='{0}/groups'.format(self.label))
        self.groups.name = self.groups.label
        try:
            self.groups.loadOldDictionary(os.path.join(path, 'dictionary.txt'), pattern=True)
        except Exception:
            logging.error('Error while reading old kinetics family dictionary from {0!r}.'.format(path))
            raise
        try:
            self.groups.loadOldTree(os.path.join(path, 'tree.txt'))
        except Exception:
            logging.error('Error while reading old kinetics family tree from {0!r}.'.format(path))
            raise

        # The old kinetics groups use rate rules (not group additivity values),
        # so we can't load the old rateLibrary.txt

        # Load the reaction recipe
        try:
            self.loadOldTemplate(os.path.join(path, 'reactionAdjList.txt'))
        except Exception:
            logging.error('Error while reading old kinetics family template/recipe from {0!r}.'.format(path))
            raise
        # Construct the forward and reverse templates
        reactants = [self.groups.entries[label] for label in self.forwardTemplate.reactants]
        if self.ownReverse:
            self.forwardTemplate = Reaction(reactants=reactants, products=reactants)
            self.reverseTemplate = None
        else:
            products = self.generateProductTemplate(reactants)
            self.forwardTemplate = Reaction(reactants=reactants, products=products)
            self.reverseTemplate = Reaction(reactants=reactants, products=products)

        self.groups.reactantNum = len(self.forwardTemplate.reactants)

        # Load forbidden structures if present
        try:
            if os.path.exists(os.path.join(path, 'forbiddenGroups.txt')):
                self.forbidden = ForbiddenStructures().loadOld(os.path.join(path, 'forbiddenGroups.txt'))
        except Exception:
            logging.error('Error while reading old kinetics family forbidden groups from {0!r}.'.format(path))
            raise

        entries = self.groups.top[:]
        for entry in self.groups.top:
            entries.extend(self.groups.descendants(entry))
        for index, entry in enumerate(entries):
            entry.index = index + 1

        self.rules = KineticsRules(label='{0}/rules'.format(self.label))
        self.rules.name = self.rules.label
        try:
            self.rules.loadOld(path, self.groups,
                               numLabels=max(len(self.forwardTemplate.reactants), len(self.groups.top)))
        except Exception:
            logging.error('Error while reading old kinetics family rules from {0!r}.'.format(path))
            raise
        self.depositories = {}

        return self

    def loadOldTemplate(self, path):
        """
        Load an old-style RMG reaction family template from the location `path`.
        """
        warnings.warn("The old kinetics databases are no longer supported and"
                      " may be removed in version 2.3.", DeprecationWarning)
        self.forwardTemplate = Reaction(reactants=[], products=[])
        self.forwardRecipe = ReactionRecipe()
        self.ownReverse = False

        ftemp = None
        # Process the template file
        try:
            ftemp = open(path, 'r')
            for line in ftemp:
                line = line.strip()
                if len(line) > 0 and line[0] == '(':
                    # This is a recipe action line
                    tokens = line.split()
                    action = [tokens[1]]
                    action.extend(tokens[2][1:-1].split(','))
                    self.forwardRecipe.addAction(action)
                elif 'thermo_consistence' in line:
                    self.ownReverse = True
                elif 'reverse' in line:
                    self.reverse = line.split(':')[1].strip()
                elif '->' in line:
                    # This is the template line
                    tokens = line.split()
                    at_arrow = False
                    for token in tokens:
                        if token == '->':
                            at_arrow = True
                        elif token != '+' and not at_arrow:
                            self.forwardTemplate.reactants.append(token)
                        elif token != '+' and at_arrow:
                            self.forwardTemplate.products.append(token)
        except IOError as e:
            logging.exception('Database template file "' + e.filename + '" not found.')
            raise
        finally:
            if ftemp: ftemp.close()

    def saveOld(self, path):
        """
        Save the old RMG kinetics groups to the given `path` on disk.
        """
        warnings.warn("The old kinetics databases are no longer supported and"
                      " may be removed in version 2.3.", DeprecationWarning)
        if not os.path.exists(path): os.mkdir(path)

        self.groups.saveOldDictionary(os.path.join(path, 'dictionary.txt'))
        self.groups.saveOldTree(os.path.join(path, 'tree.txt'))
        # The old kinetics groups use rate rules (not group additivity values),
        # so we can't save the old rateLibrary.txt
        self.saveOldTemplate(os.path.join(path, 'reactionAdjList.txt'))
        # Save forbidden structures if present
        if self.forbidden is not None:
            self.forbidden.saveOld(os.path.join(path, 'forbiddenGroups.txt'))

        self.rules.saveOld(path, self)

    def saveOldTemplate(self, path):
        """
        Save an old-style RMG reaction family template from the location `path`.
        """
        warnings.warn("The old kinetics databases are no longer supported and"
                      " may be removed in version 2.3.", DeprecationWarning)
        f_temp = open(path, 'w')

        # Write the template
        f_temp.write('{0} -> {1}\n'.format(
            ' + '.join([entry.label for entry in self.forwardTemplate.reactants]),
            ' + '.join([entry.label for entry in self.forwardTemplate.products]),
        ))
        f_temp.write('\n')

        # Write the reaction type and reverse name
        if self.ownReverse:
            f_temp.write('thermo_consistence\n')
        else:
            f_temp.write('forward\n')
            f_temp.write('reverse: {0}\n'.format(self.reverse))
        f_temp.write('\n')

        # Write the reaction recipe
        f_temp.write('Actions 1\n')
        for index, action in enumerate(self.forwardRecipe.actions):
            f_temp.write('({0}) {1:<15} {{{2}}}\n'.format(index + 1, action[0], ','.join(action[1:])))
        f_temp.write('\n')

        f_temp.close()

    def distributeTreeDistances(self):
        """
        fills in nodalDistance (the distance between an entry and its parent)
        if not already entered with the value from tree_distances associated
        with the tree the entry comes from
        """
        tree_distances = self.treeDistances
        top_labels = [i.label for i in self.groups.top]

        if len(top_labels) != len(tree_distances):
            raise ValueError('tree_distances does not have the same number of '
                             'entries as there are top nodes in the family')

        for entry_name, entry in self.groups.entries.items():
            top_entry = entry
            while not (top_entry.parent is None):  # get the top for the tree entry is in
                top_entry = top_entry.parent
            if top_entry.label in top_labels:  # filtering out product nodes
                if entry.nodalDistance is None:
                    entry.nodalDistance = tree_distances[top_entry.label]

    def load(self, path, local_context=None, global_context=None, depositoryLabels=None):
        """
        Load a kinetics database from a file located at `path` on disk.
        
        If `depositoryLabels` is a list, eg. ['training','PrIMe'], then only those
        depositories are loaded, and they are searched in that order when
        generating kinetics.
        
        If depositoryLabels is None then load 'training' first then everything else.
        If depositoryLabels is not None then load in the order specified in depositoryLabels.
        """
        local_context['recipe'] = self.loadRecipe
        local_context['template'] = self.loadTemplate
        local_context['forbidden'] = self.loadForbidden
        local_context['True'] = True
        local_context['False'] = False
        local_context['reverse'] = None
        local_context['reversible'] = None
        local_context['boundaryAtoms'] = None
        local_context['treeDistances'] = None
        local_context['reverseMap'] = None
        local_context['reactantNum'] = None
        local_context['productNum'] = None
        local_context['autoGenerated'] = False
        self.groups = KineticsGroups(label='{0}/groups'.format(self.label))
        logging.debug("Loading kinetics family groups from {0}".format(os.path.join(path, 'groups.py')))
        Database.load(self.groups, os.path.join(path, 'groups.py'), local_context, global_context)
        self.name = self.label
        self.boundaryAtoms = local_context.get('boundaryAtoms', None)
        self.treeDistances = local_context.get('treeDistances', None)
        self.reverseMap = local_context.get('reverseMap', None)

        self.reactantNum = local_context.get('reactantNum', None)
        self.productNum = local_context.get('productNum', None)

        self.autoGenerated = local_context.get('autoGenerated', False)

        if self.reactantNum:
            self.groups.reactantNum = self.reactantNum
        else:
            self.groups.reactantNum = len(self.forwardTemplate.reactants)

        # Generate the reverse template if necessary
        self.forwardTemplate.reactants = [self.groups.entries[label] for label in self.forwardTemplate.reactants]
        if self.ownReverse:
            self.forwardTemplate.products = self.forwardTemplate.reactants[:]
            self.reverseTemplate = None
            self.reverseRecipe = self.forwardRecipe.getReverse()
        else:
            self.reverse = local_context.get('reverse', None)
            self.reversible = True if local_context.get('reversible', None) is None else local_context.get('reversible', None)
            self.forwardTemplate.products = self.generateProductTemplate(self.forwardTemplate.reactants)
            if self.reversible:
                self.reverseTemplate = Reaction(reactants=self.forwardTemplate.products,
                                                products=self.forwardTemplate.reactants)
                self.reverseRecipe = self.forwardRecipe.getReverse()
                if self.reverse is None:
                    self.reverse = '{0}_reverse'.format(self.label)

        self.rules = KineticsRules(label='{0}/rules'.format(self.label))
        logging.debug("Loading kinetics family rules from {0}".format(os.path.join(path, 'rules.py')))
        self.rules.load(os.path.join(path, 'rules.py'), local_context, global_context)

        # load the groups indicated in the entry label
        for label, entries in self.rules.entries.items():
            nodes = label.split(';')
            reactants = [self.groups.entries[node] for node in nodes]
            reaction = Reaction(reactants=reactants, products=[])
            for entry in entries:
                entry.item = reaction
        self.depositories = []

        top_labels = [i.label for i in self.groups.top]
        if self.treeDistances is None:
            self.treeDistances = {top_entry: 1 for top_entry in top_labels}

        self.distributeTreeDistances()

        if depositoryLabels == 'all':
            # Load everything. This option is generally used for working with the database
            # load all the remaining depositories, in order returned by os.walk
            for root, dirs, files in os.walk(path):
                for name in dirs:
                    # if not f.endswith('.py'): continue
                    # name = f.split('.py')[0]
                    # if name not in ['groups', 'rules']:
                    f_path = os.path.join(path, name, 'reactions.py')
                    label = '{0}/{1}'.format(self.label, name)
                    depository = KineticsDepository(label=label)
                    logging.debug("Loading kinetics family depository from {0}".format(f_path))
                    depository.load(f_path, local_context, global_context)
                    self.depositories.append(depository)

            return

        if not depositoryLabels:
            # If depository labels is None or there are no depositories listed, then use the training
            # depository and add them to the RMG rate rules by default:
            depositoryLabels = ['training']
        if depositoryLabels:
            # If there are depository labels, load them in the order specified, but 
            # append the training reactions unless the user specifically declares it not
            # to be included with a '!training' flag
            if '!training' not in depositoryLabels:
                if 'training' not in depositoryLabels:
                    depositoryLabels.append('training')

        for name in depositoryLabels:
            if name == '!training':
                continue
            label = '{0}/{1}'.format(self.label, name)
            # f = name+'.py'
            f_path = os.path.join(path, name, 'reactions.py')
            if not os.path.exists(f_path):
                logging.warning("Requested depository {0} does not exist".format(f_path))
                continue
            depository = KineticsDepository(label=label)
            logging.debug("Loading kinetics family depository from {0}".format(f_path))
            depository.load(f_path, local_context, global_context)
            self.depositories.append(depository)

    def loadTemplate(self, reactants, products, ownReverse=False):
        """
        Load information about the reaction template.
        """
        self.forwardTemplate = Reaction(reactants=reactants, products=products)
        self.ownReverse = ownReverse

    def loadRecipe(self, actions):
        """
        Load information about the reaction recipe.
        """
        # Remaining lines are reaction recipe for forward reaction
        self.forwardRecipe = ReactionRecipe()
        for action in actions:
            action[0] = action[0].upper()
            valid_actions = [
                'CHANGE_BOND', 'FORM_BOND', 'BREAK_BOND', 'GAIN_RADICAL', 'LOSE_RADICAL', 'GAIN_PAIR', 'LOSE_PAIR'
            ]
            if action[0] not in valid_actions:
                raise InvalidActionError('Action {0} is not a recognized action. '
                                         'Should be one of {1}'.format(actions[0], valid_actions))
            self.forwardRecipe.addAction(action)

    def loadForbidden(self, label, group, shortDesc='', longDesc=''):
        """
        Load information about a forbidden structure.
        """
        if not self.forbidden:
            self.forbidden = ForbiddenStructures()
        self.forbidden.loadEntry(label=label, group=group, shortDesc=shortDesc, longDesc=longDesc)

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)

    def saveTrainingReactions(self, reactions, reference=None, referenceType='', shortDesc='', longDesc='', rank=3):
        """
        This function takes a list of reactions appends it to the training reactions file.  It ignores the existence of
        duplicate reactions.  
        
        The rank for each new reaction's kinetics is set to a default value of 3 unless the user specifies differently 
        for those reactions.
        
        For each entry, the long description is imported from the kinetics comment. 
        """

        if not isinstance(reference, list):
            reference = [reference] * len(reactions)
        if not isinstance(referenceType, list):
            referenceType = [referenceType] * len(reactions)
        if not isinstance(shortDesc, list):
            shortDesc = [shortDesc] * len(reactions)
        if not isinstance(longDesc, list):
            longDesc = [longDesc] * len(reactions)
        if not isinstance(rank, list):
            rank = [rank] * len(reactions)

        training_path = os.path.join(settings['database.directory'], 'kinetics', 'families',
                                     self.label, 'training')

        dictionary_path = os.path.join(training_path, 'dictionary.txt')

        # Load the old set of the species of the training reactions
        species_dict = Database().getSpecies(dictionary_path)

        # Add new unique species with labeledAtoms into species_dict
        for rxn in reactions:
            for spec in (rxn.reactants + rxn.products):
                for ex_spec in species_dict.values():
                    if ex_spec.molecule[0].getFormula() != spec.molecule[0].getFormula():
                        continue
                    else:
                        spec_labeled_atoms = spec.molecule[0].getLabeledAtoms()
                        ex_spec_labeled_atoms = ex_spec.molecule[0].getLabeledAtoms()
                        initial_map = {}
                        try:
                            for atomLabel in spec_labeled_atoms:
                                initial_map[spec_labeled_atoms[atomLabel]] = ex_spec_labeled_atoms[atomLabel]
                        except KeyError:
                            # Atom labels did not match, therefore not a match
                            continue
                        if spec.molecule[0].isIsomorphic(ex_spec.molecule[0], initial_map):
                            spec.label = ex_spec.label
                            break
                else:  # No isomorphic existing species found
                    spec_formula = spec.molecule[0].getFormula()
                    if spec_formula not in species_dict:
                        spec.label = spec_formula
                    else:
                        index = 2
                        while (spec_formula + '-{}'.format(index)) in species_dict:
                            index += 1
                        spec.label = spec_formula + '-{}'.format(index)
                    species_dict[spec.label] = spec

        training_file = open(os.path.join(training_path, 'reactions.py'), 'a')

        # get max reaction entry index from the existing training data
        try:
            depository = self.getTrainingDepository()
        except:
            logging.info('Could not find training depository in family {0}.'.format(self.label))
            logging.info('Starting a new one')
            depository = KineticsDepository()
            self.depositories.append(depository)

        if depository.entries:
            max_index = max(depository.entries.keys())
        else:
            max_index = 0

        # Add new reactions to training depository
        for i, reaction in enumerate(reactions):
            index = max_index + i + 1
            entry = Entry(
                index=index,
                label=str(reaction),
                item=reaction,
                data=reaction.kinetics,
                reference=reference[i],
                referenceType=referenceType[i],
                shortDesc=str(shortDesc[i]),
                longDesc=str(longDesc[i]),
                rank=rank[i],
            )

            # Add this entry to the loaded depository so it is immediately usable
            depository.entries[index] = entry
            # Write the entry to the reactions.py file
            self.saveEntry(training_file, entry)

        training_file.close()

        # save species to dictionary
        with open(dictionary_path, 'w') as f:
            for label in species_dict.keys():
                f.write(species_dict[label].molecule[0].toAdjacencyList(label=label, removeH=False))
                f.write('\n')

    def save(self, path):
        """
        Save the current database to the file at location `path` on disk. 
        """
        self.saveGroups(os.path.join(path, 'groups.py'))
        self.rules.save(os.path.join(path, 'rules.py'))
        for depository in self.depositories:
            self.saveDepository(depository, os.path.join(path, '{0}'.format(depository.label[len(self.label) + 1:])))

    def saveDepository(self, depository, path):
        """
        Save the given kinetics family `depository` to the location `path` on
        disk.
        """
        depository.saveDictionary(os.path.join(path, 'dictionary.txt'))
        depository.save(os.path.join(path, 'reactions.py'))

    def saveGroups(self, path):
        """
        Save the current database to the file at location `path` on disk. 
        """
        entries = self.groups.getEntriesToSave()

        # Write the header
        f = codecs.open(path, 'w', 'utf-8')
        f.write('#!/usr/bin/env python\n')
        f.write('# encoding: utf-8\n\n')
        f.write('name = "{0}/groups"\n'.format(self.name))
        f.write('shortDesc = "{0}"\n'.format(self.groups.shortDesc))
        f.write('longDesc = """\n')
        f.write(self.groups.longDesc)
        f.write('\n"""\n\n')

        # Write the template
        f.write('template(reactants=[{0}], products=[{1}], ownReverse={2})\n\n'.format(
            ', '.join(['"{0}"'.format(entry.label) for entry in self.forwardTemplate.reactants]),
            ', '.join(['"{0}"'.format(entry.label) for entry in self.forwardTemplate.products]),
            self.ownReverse))

        # Write reverse name
        if not self.ownReverse:
            if self.reverse is not None:
                f.write('reverse = "{0}"\n'.format(self.reverse))
            else:
                f.write('reverse = None\n')

        f.write('reversible = {0}\n\n'.format(self.reversible))

        if self.reverseMap is not None:
            f.write('reverseMap = {0}\n\n'.format(self.reverseMap))

        if self.reactantNum is not None:
            f.write('reactantNum = {0}\n\n'.format(self.reactantNum))
        if self.productNum is not None:
            f.write('productNum = {0}\n\n'.format(self.productNum))

        # Write the recipe
        f.write('recipe(actions=[\n')
        for action in self.forwardRecipe.actions:
            f.write('    {0!r},\n'.format(action))
        f.write('])\n\n')

        if self.boundaryAtoms:
            f.write('boundaryAtoms = ["{0}", "{1}"]'.format(self.boundaryAtoms[0], self.boundaryAtoms[1]))
            f.write('\n\n')

        # Save the entries
        for entry in entries:
            self.saveEntry(f, entry)

        # Write the tree
        if len(self.groups.top) > 0:
            f.write('tree(\n')
            f.write('"""\n')
            f.write(self.generateOldTree(self.groups.top, 1))
            f.write('"""\n')
            f.write(')\n\n')

        # Save forbidden structures, if present
        if self.forbidden is not None:
            entries = list(self.forbidden.entries.values())
            entries.sort(key=lambda x: x.label)
            for entry in entries:
                self.forbidden.saveEntry(f, entry, name='forbidden')

        f.close()

    def generateProductTemplate(self, reactants0):
        """
        Generate the product structures by applying the reaction template to
        the top-level nodes. For reactants defined by multiple structures, only
        the first is used here; it is assumed to be the most generic.
        """

        # First, generate a list of reactant structures that are actual
        # structures, rather than unions
        reactant_structures = []
        logging.log(1, "Generating template for products.")
        for reactant in reactants0:
            if isinstance(reactant, list):
                reactants = [reactant[0]]
            else:
                reactants = [reactant]

            logging.log(1, "Reactants: {0}".format(reactants))
            for s in reactants:
                logging.log(1, "Reactant {0}".format(s))
                struct = s.item
                if isinstance(struct, LogicNode):
                    all_structures = struct.getPossibleStructures(self.groups.entries)
                    logging.log(1, 'Expanding logic node {0} to {1}'.format(s, all_structures))
                    reactant_structures.append(all_structures)
                    for p in all_structures:
                        logging.log(1, p.toAdjacencyList())
                else:
                    reactant_structures.append([struct])
                    logging.log(1, struct.toAdjacencyList())

        # Second, get all possible combinations of reactant structures
        reactant_structures = getAllCombinations(reactant_structures)

        # Third, generate all possible product structures by applying the
        # recipe to each combination of reactant structures
        # Note that bimolecular products are split by labeled atoms
        product_structures = []
        for reactant_structure in reactant_structures:
            product_structure = self.applyRecipe(reactant_structure, forward=True, unique=False)
            if product_structure:
                product_structures.append(product_structure)

        # Fourth, remove duplicates from the lists
        product_structure_list = [[] for i in range(len(product_structures[0]))]
        for product_structure in product_structures:
            for i, struct in enumerate(product_structure):
                for s in product_structure_list[i]:
                    try:
                        if s.isIdentical(struct): break
                    except KeyError:
                        logging.error(struct.toAdjacencyList())
                        logging.error(s.toAdjacencyList())
                        raise
                else:
                    product_structure_list[i].append(struct)

        logging.log(1, "Unique generated product structures:")
        logging.log(1, "\n".join([p[0].toAdjacencyList() for p in product_structures]))

        # Fifth, associate structures with product template
        product_set = []
        for index, products in enumerate(product_structure_list):
            label = self.forwardTemplate.products[index]
            if len(products) == 1:
                entry = Entry(
                    label=label,
                    item=products[0],
                )
                self.groups.entries[entry.label] = entry
                product_set.append(entry)
            else:
                children = []
                counter = 0
                for product in products:
                    entry = Entry(
                        label='{0}{1:d}'.format(label, counter + 1),
                        item=product,
                    )
                    children.append(entry)
                    self.groups.entries[entry.label] = entry
                    counter += 1

                # Enter the parent of the groups as a logicOr of all the products
                entry = Entry(
                    label=label,
                    item=LogicOr([child.label for child in children], invert=False),
                    children=children,
                )
                self.groups.entries[entry.label] = entry
                # Make this entry the parent of all its children
                for child in children:
                    child.parent = entry
                counter += 1
                product_set.append(entry)

        return product_set

    def hasRateRule(self, template):
        """
        Return ``True`` if a rate rule with the given `template` currently 
        exists, or ``False`` otherwise.
        """
        return self.rules.hasRule(template)

    def getRateRule(self, template):
        """
        Return the rate rule with the given `template`. Raises a 
        :class:`ValueError` if no corresponding entry exists.
        """
        entry = self.rules.getRule(template)
        if entry is None:
            raise ValueError('No entry for template {0}.'.format(template))
        return entry

    def addKineticsRulesFromTrainingSet(self, thermoDatabase=None, trainIndices=None):
        """
        For each reaction involving real reactants and products in the training
        set, add a rate rule for that reaction.
        """
        try:
            depository = self.getTrainingDepository()
        except:
            logging.info('Could not find training depository in family {0}.'.format(self.label))
            logging.info('Must be because you turned off the training depository.')
            return

        # Determine number of parallel processes.
        from rmgpy.rmg.main import determine_procnum_from_RAM
        procnum = determine_procnum_from_RAM()

        tentries = depository.entries

        index = max([e.index for e in self.rules.getEntries()] or [0]) + 1

        entries = list(depository.entries.values())
        entries.sort(key=lambda x: x.index)

        if trainIndices is not None:
            entries = np.array(entries)
            entries = entries[trainIndices]

        reverse_entries = []
        for entry in entries:
            try:
                template = self.getReactionTemplate(entry.item)
            except UndeterminableKineticsError:
                # Some entries might be stored in the reverse direction for
                # this family; save them so we can try this
                reverse_entries.append(entry)
                continue

            tentries[entry.index].item.is_forward = True

            assert isinstance(entry.data, Arrhenius)
            data = deepcopy(entry.data)
            data.changeT0(1)

            if type(data) is Arrhenius:
                # more specific than isinstance(data,Arrhenius) because we want to exclude inherited subclasses!
                data = data.toArrheniusEP()
            elif isinstance(data, StickingCoefficient):
                data = StickingCoefficientBEP(
                    # todo: perhaps make a method StickingCoefficient.StickingCoefficientBEP
                    #  analogous to Arrhenius.toArrheniusEP
                    A=deepcopy(data.A),
                    n=deepcopy(data.n),
                    alpha=0,
                    E0=deepcopy(data.Ea),
                    Tmin=deepcopy(data.Tmin),
                    Tmax=deepcopy(data.Tmax)
                )
            elif isinstance(data, SurfaceArrhenius):
                data = SurfaceArrheniusBEP(
                    # todo: perhaps make a method SurfaceArrhenius.toSurfaceArrheniusBEP
                    #  analogous to Arrhenius.toArrheniusEP
                    A=deepcopy(data.A),
                    n=deepcopy(data.n),
                    alpha=0,
                    E0=deepcopy(data.Ea),
                    Tmin=deepcopy(data.Tmin),
                    Tmax=deepcopy(data.Tmax)
                )
            else:
                raise NotImplementedError("Unexpected training kinetics type {} for {}".format(type(data), entry))

            new_entry = Entry(
                index=index,
                label=';'.join([g.label for g in template]),
                item=Reaction(reactants=[g.item for g in template], products=[]),
                data=data,
                rank=entry.rank,
                reference=entry.reference,
                shortDesc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.shortDesc,
                longDesc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.longDesc,
            )
            new_entry.data.comment = "From training reaction {1} used for {0}".format(
                ';'.join([g.label for g in template]), entry.index)

            new_entry.data.A.value_si /= entry.item.degeneracy
            try:
                self.rules.entries[new_entry.label].append(new_entry)
            except KeyError:
                self.rules.entries[new_entry.label] = [new_entry]
            index += 1

        # Process the entries that are stored in the reverse direction of the
        # family definition
        for entry in reverse_entries:

            tentries[entry.index].item.is_forward = False

            assert isinstance(entry.data, Arrhenius)
            data = deepcopy(entry.data)
            data.changeT0(1)
            # Estimate the thermo for the reactants and products
            # trainingSet=True used later to does not allow species to match a liquid phase library
            # and get corrected thermo which will affect reverse rate calculation
            item = Reaction(reactants=[Species(molecule=[m.molecule[0].copy(deep=True)], label=m.label)
                                       for m in entry.item.reactants],
                            products=[Species(molecule=[m.molecule[0].copy(deep=True)], label=m.label)
                                      for m in entry.item.products])

            if procnum > 1:
                # If QMTP and multiprocessing write QMTP files here in parallel.
                from rmgpy.rmg.input import getInput
                quantum_mechanics = getInput('quantumMechanics')
                if quantum_mechanics:
                    quantum_mechanics.runJobs(item.reactants + item.products, procnum=procnum)

            for reactant in item.reactants:
                # Clear atom labels to avoid effects on thermo generation, ok because this is a deepcopy
                reactant.molecule[0].clearLabeledAtoms()
                reactant.generate_resonance_structures()
                reactant.thermo = thermoDatabase.getThermoData(reactant, trainingSet=True)
            for product in item.products:
                product.molecule[0].clearLabeledAtoms()
                product.generate_resonance_structures()
                product.thermo = thermoDatabase.getThermoData(product, trainingSet=True)
            # Now that we have the thermo, we can get the reverse k(T)
            item.kinetics = data
            data = item.generateReverseRateCoefficient()

            item = TemplateReaction(reactants=[m.molecule[0].copy(deep=True) for m in entry.item.products],
                                    products=[m.molecule[0].copy(deep=True) for m in entry.item.reactants])
            template = self.getReactionTemplate(item)

            item.template = self.getReactionTemplateLabels(item)
            new_degeneracy = self.calculateDegeneracy(item)

            new_entry = Entry(
                index=index,
                label=';'.join([g.label for g in template]),
                item=Reaction(reactants=[g.item for g in template],
                              products=[]),
                data=data.toArrheniusEP(),
                rank=entry.rank,
                reference=entry.reference,
                shortDesc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.shortDesc,
                longDesc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.longDesc,
            )
            new_entry.data.comment = "From training reaction {1} used for {0}".format(';'.join([g.label for g in template]), entry.index)

            new_entry.data.A.value_si /= new_degeneracy
            try:
                self.rules.entries[new_entry.label].append(new_entry)
            except KeyError:
                self.rules.entries[new_entry.label] = [new_entry]
            index += 1

    def getRootTemplate(self):
        """
        Return the root template for the reaction family. Most of the time this
        is the top-level nodes of the tree (as stored in the 
        :class:`KineticsGroups` object), but there are a few exceptions (e.g.
        R_Recombination).
        """
        if len(self.forwardTemplate.reactants) > len(self.groups.top):
            return self.forwardTemplate.reactants
        else:
            return self.groups.top

    def fillKineticsRulesByAveragingUp(self, verbose=False):
        """
        Fill in gaps in the kinetics rate rules by averaging child nodes
        recursively starting from the top level root template.
        """

        self.rules.fillRulesByAveragingUp(self.getRootTemplate(), {}, verbose)

    def applyRecipe(self, reactantStructures, forward=True, unique=True):
        """
        Apply the recipe for this reaction family to the list of
        :class:`Molecule` objects `reactantStructures`. The atoms
        of the reactant structures must already be tagged with the appropriate
        labels. Returns a list of structures corresponding to the products
        after checking that the correct number of products was produced.
        """

        # There is some hardcoding of reaction families in this function, so
        # we need the label of the reaction family for this
        label = self.label.lower()

        # Merge reactant structures into single structure
        # Also copy structures so we don't modify the originals
        # Since the tagging has already occurred, both the reactants and the
        # products will have tags
        if isinstance(reactantStructures[0], Group):
            reactant_structure = Group()
        else:
            reactant_structure = Molecule()
        for s in reactantStructures:
            reactant_structure = reactant_structure.merge(s.copy(deep=True))

        if forward:
            # Hardcoding of reaction family for peroxyl disproportionation
            # '*1' and '*2' have to be changed to '*3' and '*4' for the second reactant
            if label == 'peroxyl_disproportionation':
                identical_center_counter1 = identical_center_counter2 = 0
                for atom in reactant_structure.atoms:
                    if atom.label == '*1':
                        identical_center_counter1 += 1
                        if identical_center_counter1 > 1:
                            atom.label = '*3'
                    elif atom.label == '*2':
                        identical_center_counter2 += 1
                        if identical_center_counter2 > 1:
                            atom.label = '*4'
                msg = 'Trying to apply recipe for reaction family {}:'.format(label)
                error = False
                if identical_center_counter1 != 2:
                    msg += ' Only one occurrence of "*1" found.'
                    error = True
                if identical_center_counter2 != 2:
                    msg += ' Only one occurrence of "*2" found.'
                    error = True
                if error:
                    raise KineticsError(msg)
            # Hardcoding of reaction family for bimolecular hydroperoxide decomposition
            # '*2' has to be changed to '*4' for the second reactant and '*1' has to be
            # changed to '*6'. '*3' has to be changed to '*5' for the first reactant.
            # '*5' and '*6' do no participate in the reaction but are required for
            # relabeling in the reverse direction.
            elif label == 'bimolec_hydroperoxide_decomposition':
                identical_center_counter1 = identical_center_counter2 = identical_center_counter3 = 0
                for atom in reactant_structure.atoms:
                    if atom.label == '*1':
                        identical_center_counter1 += 1
                        if identical_center_counter1 > 1:
                            atom.label = '*6'
                    elif atom.label == '*2':
                        identical_center_counter2 += 1
                        if identical_center_counter2 > 1:
                            atom.label = '*4'
                    elif atom.label == '*3':
                        identical_center_counter3 += 1
                        if identical_center_counter3 == 1:
                            atom.label = '*5'
                msg = 'Trying to apply recipe for reaction family {}:'.format(label)
                error = False
                if identical_center_counter1 != 2:
                    msg += ' Only one occurrence of "*1" found.'
                    error = True
                if identical_center_counter2 != 2:
                    msg += ' Only one occurrence of "*2" found.'
                    error = True
                if identical_center_counter3 != 2:
                    msg += ' Only one occurrence of "*3" found.'
                    error = True
                if error:
                    raise KineticsError(msg)

            # Generate the product structure by applying the recipe
            self.forwardRecipe.applyForward(reactant_structure, unique)
        else:
            self.reverseRecipe.applyForward(reactant_structure, unique)

        # Now that we have applied the recipe, let's start calling
        # this thing the product_structure (although it's the same object in memory)
        product_structure = reactant_structure

        if not product_structure.props['validAromatic']:
            if isinstance(product_structure, Molecule):
                # For molecules, kekulize the product to redistribute bonds appropriately
                product_structure.kekulize()
            else:
                # For groups, we ignore the product template for a purely aromatic group
                # If there is an analagous aliphatic group in the family, then the product template will be identical
                # There should NOT be any families that consist solely of aromatic reactant templates
                return []

        if not forward:
            # Hardcoding of reaction family for reverse of peroxyl disproportionation
            # Labels '*3' and '*4' have to be changed back to '*1' and '*2'
            if label == 'peroxyl_disproportionation':
                for atom in product_structure.atoms:
                    if atom.label == '*3':
                        atom.label = '*1'
                    elif atom.label == '*4':
                        atom.label = '*2'
            # Hardcoding of reaction family for bimolecular hydroperoxide decomposition
            # '*5' has to be changed back to '*3', '*6' has to be changed to '*1', and
            # '*4' has to be changed to '*2'
            elif label == 'bimolec_hydroperoxide_decomposition':
                for atom in product_structure.atoms:
                    if atom.label == '*5':
                        atom.label = '*3'
                    elif atom.label == '*6':
                        atom.label = '*1'
                    elif atom.label == '*4':
                        atom.label = '*2'

        # If reaction family is its own reverse, relabel atoms
        # This allows comparison of the product species to forbidden
        #  structures which are labeled as reactants.
        # Unfortunately, this means that reaction family info is
        #  hardcoded, so this must be updated if the database changes.
        if not self.reverseTemplate:
            # Get atom labels for products
            atom_labels = {}
            for atom in product_structure.atoms:
                if atom.label != '':
                    atom_labels[atom.label] = atom

            if label == 'h_abstraction':
                # '*2' is the H that migrates
                # it moves from '*1' to '*3'
                atom_labels['*1'].label = '*3'
                atom_labels['*3'].label = '*1'

            elif label == 'intra_h_migration':
                # '*3' is the H that migrates
                # swap the two ends between which the H moves
                atom_labels['*1'].label = '*2'
                atom_labels['*2'].label = '*1'
                # reverse all the atoms in the chain between *1 and *2
                highest = len(atom_labels)
                if highest > 4:
                    # swap *4 with *5
                    atom_labels['*4'].label = '*5'
                    atom_labels['*5'].label = '*4'
                if highest > 6:
                    # swap *6 with the highest, etc.
                    for i in range(6, highest + 1):
                        atom_labels['*{0:d}'.format(i)].label = '*{0:d}'.format(6 + highest - i)

            elif label == 'intra_ene_reaction':
                # Labels for nodes are swapped
                atom_labels['*1'].label = '*2'
                atom_labels['*2'].label = '*1'
                atom_labels['*3'].label = '*5'
                atom_labels['*5'].label = '*3'

            elif label == '6_membered_central_c-c_shift':
                # Labels for nodes are swapped
                atom_labels['*1'].label = '*3'
                atom_labels['*3'].label = '*1'
                atom_labels['*4'].label = '*6'
                atom_labels['*6'].label = '*4'

            elif label == '1,2_shiftc':
                # Labels for nodes are swapped
                atom_labels['*2'].label = '*3'
                atom_labels['*3'].label = '*2'

            elif label == 'intra_r_add_exo_scission':
                # Labels for nodes are swapped
                atom_labels['*1'].label = '*3'
                atom_labels['*3'].label = '*1'

            elif label == 'intra_substitutions_isomerization':
                # Swap *2 and *3
                atom_labels['*2'].label = '*3'
                atom_labels['*3'].label = '*2'

        if not forward:
            template = self.reverseTemplate
            product_num = self.reactantNum or len(template.products)
        else:
            template = self.forwardTemplate
            product_num = self.productNum or len(template.products)

        # Split product structure into multiple species if necessary
        product_structures = product_structure.split()

        # Make sure we've made the expected number of products
        if product_num != len(product_structures):
            # We have a different number of products than expected by the template.
            # By definition this means that the template is not a match, so
            # we return None to indicate that we could not generate the product
            # structures
            # We need to think this way in order to distinguish between
            # intermolecular and intramolecular versions of reaction families,
            # which will have very different kinetics
            # Unfortunately this may also squash actual errors with malformed
            # reaction templates
            return None

        # Remove vdW bonds
        for struct in product_structures:
            struct.removeVanDerWaalsBonds()

        # Make sure we don't create a different net charge between reactants and products
        reactant_net_charge = product_net_charge = 0
        for struc in reactantStructures:
            struc.update()
            reactant_net_charge += struc.getNetCharge()

        for struct in product_structures:
            # If product structures are Molecule objects, update their atom types
            # If product structures are Group objects and the reaction is in certain families
            # (families with charged substances), the charge of structures will be updated
            if isinstance(struct, Molecule):
                struct.update()
            elif isinstance(struct, Group):
                struct.resetRingMembership()
                if label in ['1,2_insertion_co', 'r_addition_com', 'co_disproportionation',
                             'intra_no2_ono_conversion', 'lone_electron_pair_bond',
                             '1,2_nh3_elimination', '1,3_nh3_elimination']:
                    struct.update_charge()
            else:
                raise TypeError('Expecting Molecule or Group object, not {0}'.format(struct.__class__.__name__))
            product_net_charge += struc.getNetCharge()
        if reactant_net_charge != product_net_charge:
            logging.debug(
                'The net charge of the reactants {0} differs from the net charge of the products {1} in reaction '
                'family {2}. Not generating this reaction.'.format(reactant_net_charge, product_net_charge, self.label))
            return None
        # The following check should be removed once RMG can process charged species
        # This is applied only for :class:Molecule (not for :class:Group which is allowed to have a nonzero net charge)
        if any([structure.getNetCharge() for structure in reactantStructures + product_structures]) \
                and isinstance(struc, Molecule):
            logging.debug(
                'A net charged species was formed when reacting {0} to form {1} in reaction family {2}. Not generating '
                'this reaction.'.format(reactant_net_charge, product_net_charge, self.label))
            return None

        # If there are two product structures, place the one containing '*1' first
        if len(product_structures) == 2:
            if not product_structures[0].containsLabeledAtom('*1') and \
                    product_structures[1].containsLabeledAtom('*1'):
                product_structures.reverse()
        # If there are three product structures, sort them based on the lowest number label in each structure
        elif len(product_structures) == 3:
            lowest_labels = []
            for struct in product_structures:
                # Extract digits from labels and convert others (e.g., "*") to empty strings
                labels = [''.join(c for c in label if c.isdigit()) for label in struct.getLabeledAtoms().keys()]
                # Convert digits to integers and remove empty strings
                labels = [int(label) for label in labels if label]
                lowest_labels.append(min(labels))
            product_structures = [s for _, s in sorted(zip(lowest_labels, product_structures))]

        # Return the product structures
        return product_structures

    def __generateProductStructures(self, reactantStructures, maps, forward):
        """
        For a given set of `reactantStructures` and a given set of `maps`,
        generate and return the corresponding product structures. The
        `reactantStructures` parameter should be given in the order the
        reactants are stored in the reaction family template. The `maps`
        parameter is a list of mappings of the top-level tree node of each
        *template* reactant to the corresponding *structure*. This function
        returns a list of the product structures.
        """

        # Clear any previous atom labeling from all reactant structures
        for struct in reactantStructures:
            struct.clearLabeledAtoms()

        # Tag atoms with labels
        for m in maps:
            for reactant_atom, template_atom in m.items():
                reactant_atom.label = template_atom.label

        # Check that reactant structures are allowed in this family
        # If not, then stop
        for struct in reactantStructures:
            if self.isMoleculeForbidden(struct):
                raise ForbiddenStructureException()

        # Generate the product structures by applying the forward reaction recipe
        try:
            product_structures = self.applyRecipe(reactantStructures, forward=forward)
            if not product_structures:
                return None
        except (InvalidActionError, KekulizationError):
            # If unable to apply the reaction recipe, then return no product structures
            return None
        except ActionError:
            logging.error('Could not generate product structures for reaction family {0} in {1} '
                          'direction'.format(self.label, 'forward' if forward else 'reverse'))
            logging.info('Reactant structures:')
            for struct in reactantStructures:
                logging.info('{0}\n{1}\n'.format(struct, struct.toAdjacencyList()))
            raise

        # Apply the generated species constraints (if given)
        for struct in product_structures:
            if self.isMoleculeForbidden(struct):
                raise ForbiddenStructureException()
            if failsSpeciesConstraints(struct):
                raise ForbiddenStructureException()

        return product_structures

    def isMoleculeForbidden(self, molecule):
        """
        Return ``True`` if the molecule is forbidden in this family, or
        ``False`` otherwise. 
        """

        # check family-specific forbidden structures 
        if self.forbidden is not None and self.forbidden.isMoleculeForbidden(molecule):
            return True

        return False

    def __createReaction(self, reactants, products, is_forward):
        """
        Create and return a new :class:`Reaction` object containing the
        provided `reactants` and `products` as lists of :class:`Molecule`
        objects.
        """

        # Make sure the products are in fact different than the reactants
        if same_species_lists(reactants, products):
            return None

        # Create and return template reaction object
        reaction = TemplateReaction(
            reactants=reactants if is_forward else products,
            products=products if is_forward else reactants,
            degeneracy=1,
            reversible=self.reversible,
            family=self.label,
            is_forward=is_forward,
        )

        # Store the labeled atoms so we can recover them later
        # (e.g. for generating reaction pairs and templates)
        labeled_atoms = []
        for reactant in reaction.reactants:
            for label, atom in reactant.getLabeledAtoms().items():
                labeled_atoms.append((label, atom))
        reaction.labeledAtoms = labeled_atoms

        return reaction

    def __matchReactantToTemplate(self, reactant, templateReactant):
        """
        Return a complete list of the mappings if the provided reactant 
        matches the provided template reactant, or an empty list if not.
        """

        if isinstance(templateReactant, list):
            templateReactant = templateReactant[0]
        if isinstance(templateReactant, Entry):
            struct = templateReactant.item
        else:
            struct = templateReactant

        reactant_contains_surface_site = reactant.containsSurfaceSite()

        if isinstance(struct, LogicNode):
            mappings = []
            for child_structure in struct.getPossibleStructures(self.groups.entries):
                if child_structure.containsSurfaceSite() != reactant_contains_surface_site:
                    # An adsorbed template can't match a gas-phase species and vice versa
                    continue
                mappings.extend(reactant.findSubgraphIsomorphisms(child_structure))
            return mappings
        elif isinstance(struct, Group):
            if struct.containsSurfaceSite() != reactant_contains_surface_site:
                # An adsorbed template can't match a gas-phase species and vice versa
                return []
            return reactant.findSubgraphIsomorphisms(struct)
        else:
            raise NotImplementedError("Not expecting template of type {}".format(type(struct)))

    def generateReactions(self, reactants, products=None, prod_resonance=True):
        """
        Generate all reactions between the provided list of one, two, or three
        `reactants`, which should be either single :class:`Molecule` objects
        or lists of same. Does not estimate the kinetics of these reactions
        at this time. Returns a list of :class:`TemplateReaction` objects
        using :class:`Molecule` objects for both reactants and products
        The reactions are constructed such that the forward direction is
        consistent with the template of this reaction family.

        Args:
            reactants (list):                List of Molecules to react.
            products (list, optional):       List of Molecules or Species of desired product structures.
            prod_resonance (bool, optional): Flag to generate resonance structures for product checking.
                Defaults to True, resonance structures are compared.

        Returns:
            List of all reactions containing Molecule objects with the
            specified reactants and products within this family.
            Degenerate reactions are returned as separate reactions.
        """
        reaction_list = []

        # Forward direction (the direction in which kinetics is defined)
        reaction_list.extend(
            self.__generateReactions(reactants, products=products, forward=True, prod_resonance=prod_resonance))

        if not self.ownReverse and self.reversible:
            # Reverse direction (the direction in which kinetics is not defined)
            reaction_list.extend(
                self.__generateReactions(reactants, products=products, forward=False, prod_resonance=prod_resonance))

        return reaction_list

    def addReverseAttribute(self, rxn, react_non_reactive=True):
        """
        For rxn (with species' objects) from families with ownReverse, this method adds a `reverse`
        attribute that contains the reverse reaction information (like degeneracy)

        Returns `True` if successful and `False` if the reverse reaction is forbidden.
        Will raise a `KineticsError` if unsuccessful for other reasons.
        """
        if self.ownReverse and all([spc.has_reactive_molecule() for spc in rxn.products]):
            # Check if the reactants are the same
            same_reactants = 0
            if len(rxn.products) == 2 and rxn.products[0].isIsomorphic(rxn.products[1]):
                same_reactants = 2
            elif len(rxn.products) == 3:
                same_01 = rxn.products[0].isIsomorphic(rxn.products[1])
                same_02 = rxn.products[0].isIsomorphic(rxn.products[2])
                if same_01 and same_02:
                    same_reactants = 3
                elif same_01 or same_02:
                    same_reactants = 2
                elif rxn.products[1].isIsomorphic(rxn.products[2]):
                    same_reactants = 2

            ensure_independent_atom_ids(rxn.products)

            reaction_list = self.__generateReactions([spc.molecule for spc in rxn.products],
                                                     products=rxn.reactants, forward=True,
                                                     react_non_reactive=react_non_reactive)
            reactions = find_degenerate_reactions(reaction_list, same_reactants, kinetics_family=self)
            if len(reactions) == 0:
                logging.error("Expecting one matching reverse reaction, not zero in reaction family {0} for "
                              "forward reaction {1}.\n".format(self.label, str(rxn)))
                logging.error("There is likely a bug in the RMG-database kinetics reaction family involving a "
                              "missing group, missing atomlabels, forbidden groups, etc.")
                for reactant in rxn.reactants:
                    logging.info("Reactant")
                    logging.info(reactant.toAdjacencyList())
                for product in rxn.products:
                    logging.info("Product")
                    logging.info(product.toAdjacencyList())
                logging.error("Debugging why no reaction was found...")
                logging.error("Checking whether the family's forbidden species have affected reaction generation...")
                # Set family's forbidden structures to empty for now to see if reaction gets generated...
                # Note that it is not necessary to check global forbidden structures, because this reaction
                # would not have been formed in the first place.
                temp_object = self.forbidden
                self.forbidden = ForbiddenStructures()  # Initialize with empty one
                try:
                    reaction_list = self.__generateReactions([spc.molecule for spc in rxn.products],
                                                            products=rxn.reactants, forward=True,
                                                            react_non_reactive=react_non_reactive)
                    reactions = find_degenerate_reactions(reaction_list, same_reactants, kinetics_family=self)
                finally:
                    self.forbidden = temp_object
                if (len(reactions) == 1 or
                        (len(reactions) > 1 and
                         all([reactions[0].isIsomorphic(other, checkTemplateRxnProducts=True)
                              for other in reactions]))):
                    logging.error("Error was fixed, the product is a forbidden structure when used as a reactant "
                                  "in the reverse direction.")
                    # This reaction should be forbidden in the forward direction as well
                    return False
                else:
                    logging.error("Still experiencing error: Expecting one matching reverse reaction, not {0} in "
                                  "reaction family {1} for forward reaction {2}."
                                  "\n".format(len(reactions), self.label, str(rxn)))
                    raise KineticsError("Did not find reverse reaction in reaction family {0} for reaction "
                                        "{1}.".format(self.label, str(rxn)))
            elif (len(reactions) > 1 and
                    not all([reactions[0].isIsomorphic(other, strict=False, checkTemplateRxnProducts=True)
                             for other in reactions])):
                logging.error("Expecting one matching reverse reaction. Recieved {0} reactions with "
                              "multiple non-isomorphic ones in reaction family {1} for "
                              "forward reaction {2}.\n".format(len(reactions), self.label, str(rxn)))
                logging.info("Found the following reverse reactions")
                for rxn0 in reactions:
                    logging.info(str(rxn0))
                    for reactant in rxn0.reactants:
                        logging.info("Reactant")
                        logging.info(reactant.toAdjacencyList())
                    for product in rxn0.products:
                        logging.info("Product")
                        logging.info(product.toAdjacencyList())
                raise KineticsError("Found multiple reverse reactions in reaction family {0} for reaction {1}, likely "
                                    "due to inconsistent resonance structure generation".format(self.label, str(rxn)))
            else:
                rxn.reverse = reactions[0]
                return True

    def calculateDegeneracy(self, reaction):
        """
        For a `reaction`  with `Molecule` or `Species` objects given in the direction in which
        the kinetics are defined, compute the reaction-path degeneracy.

        This method by default adjusts for double counting of identical reactants. 
        This should only be adjusted once per reaction. To not adjust for 
        identical reactants (since you will be reducing them later in the algorithm), add
        `ignoreSameReactants= True` to this method.
        """
        # Check if the reactants are the same
        # If they refer to the same memory address, then make a deep copy so
        # they can be manipulated independently
        reactants = reaction.reactants
        same_reactants = 0
        if len(reactants) == 2:
            if reactants[0] is reactants[1]:
                reactants[1] = reactants[1].copy(deep=True)
                same_reactants = 2
            elif reactants[0].isIsomorphic(reactants[1]):
                same_reactants = 2
        elif len(reactants) == 3:
            same_01 = reactants[0] is reactants[1]
            same_02 = reactants[0] is reactants[2]
            if same_01 and same_02:
                same_reactants = 3
                reactants[1] = reactants[1].copy(deep=True)
                reactants[2] = reactants[2].copy(deep=True)
            elif same_01:
                same_reactants = 2
                reactants[1] = reactants[1].copy(deep=True)
            elif same_02:
                same_reactants = 2
                reactants[2] = reactants[2].copy(deep=True)
            elif reactants[1] is reactants[2]:
                same_reactants = 2
                reactants[2] = reactants[2].copy(deep=True)
            else:
                same_01 = reactants[0].isIsomorphic(reactants[1])
                same_02 = reactants[0].isIsomorphic(reactants[2])
                if same_01 and same_02:
                    same_reactants = 3
                elif same_01 or same_02:
                    same_reactants = 2
                elif reactants[1].isIsomorphic(reactants[2]):
                    same_reactants = 2

        # Label reactant atoms for proper degeneracy calculation
        ensure_independent_atom_ids(reactants, resonance=True)
        molecule_combos = generate_molecule_combos(reactants)

        reactions = []
        for combo in molecule_combos:
            reactions.extend(self.__generateReactions(combo, products=reaction.products, forward=True,
                                                      react_non_reactive=True))

        # remove degenerate reactions
        reactions = find_degenerate_reactions(reactions, same_reactants, template=reaction.template,
                                              kinetics_family=self)

        # log issues
        if len(reactions) != 1:
            for reactant in reaction.reactants:
                logging.error("Reactant: {0!r}".format(reactant))
            for product in reaction.products:
                logging.error("Product: {0!r}".format(product))
            raise KineticsError(('Unable to calculate degeneracy for reaction {0} '
                                 'in reaction family {1}. Expected 1 reaction '
                                 'but generated {2}').format(reaction, self.label, len(reactions)))
        return reactions[0].degeneracy

    def __generateReactions(self, reactants, products=None, forward=True, prod_resonance=True,
                            react_non_reactive=False):
        """
        Generate a list of all the possible reactions of this family between
        the list of `reactants`. The number of reactants provided must match
        the number of reactants expected by the template, or this function
        will return an empty list. Each item in the list of reactants should
        be a list of :class:`Molecule` objects, each representing a resonance
        structure of the species of interest.
        
        This method returns all reactions, and degenerate reactions can then be
        found using `rmgpy.data.kinetics.common.find_degenerate_reactions`.

        Args:
            reactants:          List of Molecules to react
            products:           List of Molecules or Species of desired product structures (optional)
            forward:            Flag to indicate whether the forward or reverse template should be applied (optional)
                                Default is True, forward template is used
            prod_resonance:     Flag to generate resonance structures for product checking (optional)
                                Default is True, resonance structures are compared
            react_non_reactive: Flag to generate reactions between unreactive molecules (optional)
                                Default is False, reactions involving unreactive molecules are not generated

        Returns:
            List of all reactions containing Molecule objects with the
                specified reactants and products within this family.
            Degenerate reactions are returned as separate reactions.
        """

        rxn_list = []

        # Wrap each reactant in a list if not already done (this is done to 
        # allow for passing multiple resonance structures for each molecule)
        # This also makes a copy of the reactants list so we don't modify the
        # original
        reactants = [reactant if isinstance(reactant, list) else [reactant] for reactant in reactants]

        if forward:
            template = self.forwardTemplate
            reactant_num = self.reactantNum
        elif self.reverseTemplate is None:
            return []
        else:
            template = self.reverseTemplate
            reactant_num = self.productNum

        if self.autoGenerated and reactant_num != len(reactants):
            return []

        if len(reactants) > len(template.reactants):
            # if the family has one template and is bimolecular split template into multiple reactants
            try:
                grps = template.reactants[0].item.split()
                template_reactants = []
                for grp in grps:
                    template_reactants.append(grp)
            except AttributeError:
                template_reactants = [x.item for x in template.reactants]
        else:
            template_reactants = [x.item for x in template.reactants]

        # Unimolecular reactants: A --> products
        if len(reactants) == 1 and len(template_reactants) == 1:

            # Iterate over all resonance isomers of the reactant
            for molecule in reactants[0]:
                if molecule.reactive or react_non_reactive:  # don't react non representative resonance isomers unless
                    # explicitly desired (e.g., when called from calculateDegeneracy)
                    mappings = self.__matchReactantToTemplate(molecule, template_reactants[0])
                    for mapping in mappings:
                        reactant_structures = [molecule]
                        try:
                            product_structures = self.__generateProductStructures(reactant_structures,
                                                                                  [mapping], forward)
                        except ForbiddenStructureException:
                            pass
                        else:
                            if product_structures is not None:
                                rxn = self.__createReaction(reactant_structures, product_structures, forward)
                                if rxn: rxn_list.append(rxn)

        # Bimolecular reactants: A + B --> products
        elif len(reactants) == 2 and len(template_reactants) == 2:

            molecules_a = reactants[0]
            molecules_b = reactants[1]

            # ToDo: try to remove this hard-coding of reaction family name..
            if 'adsorption' in self.label.lower() and forward:
                if molecules_a[0].containsSurfaceSite() and molecules_b[0].containsSurfaceSite():
                    # Can't adsorb something that's already adsorbed.
                    # Both reactants either contain or are a surface site.
                    return []

            # Iterate over all resonance isomers of the reactant
            for molecule_a in molecules_a:
                for molecule_b in molecules_b:
                    if (molecule_a.reactive and molecule_b.reactive) or react_non_reactive:

                        # Reactants stored as A + B
                        mappings_a = self.__matchReactantToTemplate(molecule_a, template_reactants[0])
                        mappings_b = self.__matchReactantToTemplate(molecule_b, template_reactants[1])

                        # Iterate over each pair of matches (A, B)
                        for map_a in mappings_a:
                            for map_b in mappings_b:
                                # Reverse the order of reactants in case we have a family with only one reactant tree
                                # that can produce different products depending on the order of reactants
                                reactant_structures = [molecule_b, molecule_a]
                                try:
                                    product_structures = self.__generateProductStructures(reactant_structures,
                                                                                          [map_b, map_a], forward)
                                except ForbiddenStructureException:
                                    pass
                                else:
                                    if product_structures is not None:
                                        rxn = self.__createReaction(reactant_structures, product_structures, forward)
                                        if rxn: rxn_list.append(rxn)

                        # Only check for swapped reactants if they are different
                        if reactants[0] is not reactants[1]:

                            # Reactants stored as B + A
                            mappings_a = self.__matchReactantToTemplate(molecule_a, template_reactants[1])
                            mappings_b = self.__matchReactantToTemplate(molecule_b, template_reactants[0])

                            # Iterate over each pair of matches (A, B)
                            for map_a in mappings_a:
                                for map_b in mappings_b:
                                    reactant_structures = [molecule_a, molecule_b]
                                    try:
                                        product_structures = self.__generateProductStructures(reactant_structures,
                                                                                              [map_a, map_b], forward)
                                    except ForbiddenStructureException:
                                        pass
                                    else:
                                        if product_structures is not None:
                                            rxn = self.__createReaction(reactant_structures, product_structures,
                                                                        forward)
                                            if rxn:
                                                rxn_list.append(rxn)

        # Termolecular reactants: A + B + C --> products
        elif len(reactants) == 2 and len(template_reactants) == 3:
            """
            Two reactants but a termolecular template.
            Could be A + X + X <=> BX + CX (dissociative adsorption)
            or A + X + X <=> AXX (bidentate adsorption)
            in which case, if one of the two reactants is an X
            then we have a match and can just use it twice.
            """
            template_sites = [r for r in template_reactants if r.isSurfaceSite()]
            if len(template_sites) == 2:
                # Two surface sites in template. If there's a site in the reactants, use it twice.
                if reactants[0][0].isSurfaceSite() and not reactants[1][0].isSurfaceSite():
                    site1 = reactants[0][0]
                    site2 = deepcopy(reactants[0][0])
                    adsorbate_molecules = reactants[1]
                    reactants.append([site2])
                elif reactants[1][0].isSurfaceSite() and not reactants[0][0].isSurfaceSite():
                    site1 = reactants[1][0]
                    site2 = deepcopy(reactants[1][0])
                    adsorbate_molecules = reactants[0]
                    reactants.append([site2])
                else:
                    # No reaction with these reactants in this template
                    return []

                if adsorbate_molecules[0].containsSurfaceSite():
                    # An adsorbed molecule can't adsorb again
                    return []

                for r in template_reactants:
                    if not r.isSurfaceSite():
                        template_adsorbate = r
                        break
                else:
                    raise KineticsError("Couldn't find non-site in template {0!r}".format(template))

                mappings_a = self.__matchReactantToTemplate(site1, template_sites[0])
                mappings_b = self.__matchReactantToTemplate(site2, template_sites[1])
                for adsorbateMolecule in adsorbate_molecules:
                    mappings_c = self.__matchReactantToTemplate(adsorbateMolecule, template_adsorbate)
                    for map_a, map_b, map_c in itertools.product(mappings_a, mappings_b, mappings_c):
                        reactant_structures = [site1, site2, adsorbateMolecule]
                        # should be in same order as reaction template recipe?
                        try:
                            product_structures = self.__generateProductStructures(reactant_structures,
                                                                                  [map_a, map_b, map_c],
                                                                                  forward)
                        except ForbiddenStructureException:
                            pass
                        else:
                            if product_structures is not None:
                                rxn = self.__createReaction(reactant_structures, product_structures, forward)
                                if rxn:
                                    rxn_list.append(rxn)
            else:
                # __generateReactions was called with mismatched number of reactants and templates
                return []

        elif len(reactants) == 3 and len(template_reactants) == 3:
            """
            This could be a surface reaction
                A + X + X <=> BX + CX  (dissociative adsorption)
                A + X + X <=> AXX      (bidentate adsorption)
            or a termolecular gas phase reaction
                A + B + C <=> stuff
            We check the two scenarios in that order.
            """
            template_sites = [r for r in template_reactants if r.isSurfaceSite()]
            if len(template_sites) == 2:
                """
                Three reactants and a termolecular template.
                Could be A + X + X <=> BX + CX (dissociative adsorption)
                or A + X + X <=> AXX (bidentate adsorption)
                that was first found in the reverse direction
                and so is being passed in with all three reactants identified.
                """
                # Should be 2 surface sites in reactants too.
                # Find them, and find mappings of the other
                m1, m2, m3 = (r[0] for r in reactants)
                if m1.isSurfaceSite() and m2.isSurfaceSite() and not m3.isSurfaceSite():
                    site1, site2 = m1, m2
                    adsorbate_molecules = reactants[2]
                elif m1.isSurfaceSite() and not m2.isSurfaceSite() and m3.isSurfaceSite():
                    site1, site2 = m1, m3
                    adsorbate_molecules = reactants[1]
                elif not m1.isSurfaceSite() and m2.isSurfaceSite() and m3.isSurfaceSite():
                    site1, site2 = m2, m3
                    adsorbate_molecules = reactants[0]
                else:
                    # Three reactants not containing two surface sites
                    return []

                if adsorbate_molecules[0].containsSurfaceSite():
                    # An adsorbed molecule can't adsorb again
                    return []

                for r in template_reactants:
                    if not r.isSurfaceSite():
                        template_adsorbate = r
                        break
                else:
                    raise KineticsError("Couldn't find non-site in template {0!r}".format(template))

                mappings_a = self.__matchReactantToTemplate(site1, template_sites[0])
                mappings_b = self.__matchReactantToTemplate(site2, template_sites[1])
                for adsorbateMolecule in adsorbate_molecules:
                    mappings_c = self.__matchReactantToTemplate(adsorbateMolecule, template_adsorbate)
                    # this just copied/pasted from above - not checked
                    for map_a, map_b, map_c in itertools.product(mappings_a, mappings_b, mappings_c):
                        reactant_structures = [site1, site2, adsorbateMolecule]
                        try:
                            product_structures = self.__generateProductStructures(reactant_structures,
                                                                                  [map_a, map_b, map_c],
                                                                                  forward)
                        except ForbiddenStructureException:
                            pass
                        else:
                            if product_structures is not None:
                                rxn = self.__createReaction(reactant_structures, product_structures, forward)
                                if rxn:
                                    rxn_list.append(rxn)

            else:
                """
                Not a bidentate surface reaction, just a gas-phase
                Trimolecular reactants: A + B + C --> products
                """
                molecules_a = reactants[0]
                molecules_b = reactants[1]
                molecules_c = reactants[2]

                # Iterate over all resonance isomers of the reactants
                for molecule_a in molecules_a:
                    for molecule_b in molecules_b:
                        for molecule_c in molecules_c:

                            def generate_products_and_reactions(order):
                                """
                                order = (0, 1, 2) corresponds to reactants stored as A + B + C, etc.
                                """
                                _mappings_a = self.__matchReactantToTemplate(molecule_a, template_reactants[order[0]])
                                _mappings_b = self.__matchReactantToTemplate(molecule_b, template_reactants[order[1]])
                                _mappings_c = self.__matchReactantToTemplate(molecule_c, template_reactants[order[2]])

                                # Iterate over each pair of matches (A, B, C)
                                for _map_a in _mappings_a:
                                    for _map_b in _mappings_b:
                                        for _map_c in _mappings_c:
                                            _reactantStructures = [molecule_a, molecule_b, molecule_c]
                                            _maps = [_map_a, _map_b, _map_c]
                                            # Reorder reactants in case we have a family with fewer reactant trees than
                                            # reactants and different reactant orders can produce different products
                                            _reactantStructures = [_reactantStructures[_i] for _i in order]
                                            _maps = [_maps[_i] for _i in order]
                                            try:
                                                _productStructures = self.__generateProductStructures(
                                                    _reactantStructures,
                                                    _maps,
                                                    forward)
                                            except ForbiddenStructureException:
                                                pass
                                            else:
                                                if _productStructures is not None:
                                                    _rxn = self.__createReaction(_reactantStructures,
                                                                                 _productStructures,
                                                                                 forward)
                                                    if _rxn:
                                                        rxn_list.append(_rxn)

                            # Reactants stored as A + B + C
                            generate_products_and_reactions((0, 1, 2))

                            # Only check for swapped reactants if they are different
                            if reactants[1] is not reactants[2]:
                                # Reactants stored as A + C + B
                                generate_products_and_reactions((0, 2, 1))
                            if reactants[0] is not reactants[1]:
                                # Reactants stored as B + A + C
                                generate_products_and_reactions((1, 0, 2))
                            if reactants[0] is not reactants[2]:
                                # Reactants stored as C + B + A
                                generate_products_and_reactions((2, 1, 0))
                                if reactants[0] is not reactants[1] and reactants[1] is not reactants[2]:
                                    # Reactants stored as C + A + B
                                    generate_products_and_reactions((2, 0, 1))
                                    # Reactants stored as B + C + A
                                    generate_products_and_reactions((1, 2, 0))

        # ToDo: try to remove this hard-coding of reaction family name..
        if not forward and 'adsorption' in self.label.lower():
            # Desorption should have desorbed something (else it was probably bidentate)
            # so delete reactions that don't make a gas-phase desorbed product
            pruned_list = []
            for reaction in rxn_list:
                for reactant in reaction.reactants:
                    if not reactant.containsSurfaceSite():
                        # found a desorbed species, we're ok
                        pruned_list.append(reaction)
                        break
                else:  # didn't break, so all species still adsorbed
                    logging.debug("Removing {0} reaction {1!s} with no desorbed species".format(self.label, reaction))
                    continue  # to next reaction immediately
            rxn_list = pruned_list

        # If products is given, remove reactions from the reaction list that
        # don't generate the given products
        if products is not None:
            rxn_list0 = rxn_list[:]
            rxn_list = []
            for reaction in rxn_list0:
                products0 = reaction.products if forward else reaction.reactants
                # Only keep reactions which give the requested products
                # If prod_resonance=True, then use strict=False to consider all resonance structures
                if same_species_lists(products, products0, strict=not prod_resonance):
                    rxn_list.append(reaction)

        # Determine the reactant-product pairs to use for flux analysis
        # Also store the reaction template (useful so we can easily get the kinetics later)
        for reaction in rxn_list:

            # Restore the labeled atoms long enough to generate some metadata
            for reactant in reaction.reactants:
                reactant.clearLabeledAtoms()
            for label, atom in reaction.labeledAtoms:
                if isinstance(atom, list):
                    for atm in atom:
                        atm.label = label
                else:
                    atom.label = label

            # Generate metadata about the reaction that we will need later
            reaction.pairs = self.getReactionPairs(reaction)
            reaction.template = self.getReactionTemplateLabels(reaction)

            # Unlabel the atoms for both reactants and products
            for species in itertools.chain(reaction.reactants, reaction.products):
                species.clearLabeledAtoms()

            # We're done with the labeled atoms, so delete the attribute
            del reaction.labeledAtoms

            # Mark reaction reversibility
            reaction.reversible = self.reversible

        # This reaction list has only checked for duplicates within itself, not
        # with the global list of reactions
        return rxn_list

    def getReactionPairs(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, return the reactant-product pairs to use when
        performing flux analysis.
        """
        pairs = []
        if len(reaction.reactants) == 1 or len(reaction.products) == 1:
            # When there is only one reactant (or one product), it is paired 
            # with each of the products (reactants)
            for reactant in reaction.reactants:
                for product in reaction.products:
                    pairs.append([reactant, product])
        elif self.label.lower() == 'h_abstraction':
            # Hardcoding for hydrogen abstraction: pair the reactant containing
            # *1 with the product containing *3 and vice versa
            assert len(reaction.reactants) == len(reaction.products) == 2
            if reaction.reactants[0].containsLabeledAtom('*1'):
                if reaction.products[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
            elif reaction.reactants[1].containsLabeledAtom('*1'):
                if reaction.products[1].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
        elif self.label.lower() in ['disproportionation', 'co_disproportionation', 'korcek_step1_cat']:
            # Hardcoding for disproportionation, co_disproportionation, korcek_step1_cat:
            # pair the reactant containing *1 with the product containing *1
            assert len(reaction.reactants) == len(reaction.products) == 2
            if reaction.reactants[0].containsLabeledAtom('*1'):
                if reaction.products[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
            elif reaction.reactants[1].containsLabeledAtom('*1'):
                if reaction.products[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
        elif self.label.lower() in ['substitution_o', 'substitutions']:
            # Hardcoding for Substitution_O: pair the reactant containing
            # *2 with the product containing *3 and vice versa
            assert len(reaction.reactants) == len(reaction.products) == 2
            if reaction.reactants[0].containsLabeledAtom('*2'):
                if reaction.products[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
            elif reaction.reactants[1].containsLabeledAtom('*2'):
                if reaction.products[1].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
        elif self.label.lower() == 'baeyer-villiger_step1_cat':
            # Hardcoding for Baeyer-Villiger_step1_cat: pair the two reactants
            # with the Criegee intermediate and pair the catalyst with itself
            assert len(reaction.reactants) == 3 and len(reaction.products) == 2
            if reaction.reactants[0].containsLabeledAtom('*5'):
                if reaction.products[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                    pairs.append([reaction.reactants[2], reaction.products[0]])
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                    pairs.append([reaction.reactants[2], reaction.products[1]])
                    pairs.append([reaction.reactants[0], reaction.products[0]])
            elif reaction.reactants[1].containsLabeledAtom('*5'):
                if reaction.products[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[2], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[2], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
            elif reaction.reactants[2].containsLabeledAtom('*5'):
                if reaction.products[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                    pairs.append([reaction.reactants[2], reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                    pairs.append([reaction.reactants[2], reaction.products[0]])
        elif self.label.lower() == 'baeyer-villiger_step2_cat':
            # Hardcoding for Baeyer-Villiger_step2_cat: pair the Criegee
            # intermediate with the two products and the catalyst with itself
            assert len(reaction.reactants) == 2 and len(reaction.products) == 3
            if reaction.products[0].containsLabeledAtom('*7'):
                if reaction.reactants[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[0], reaction.products[2]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                elif reaction.reactants[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[2]])
                    pairs.append([reaction.reactants[0], reaction.products[0]])
            elif reaction.products[1].containsLabeledAtom('*7'):
                if reaction.reactants[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[0], reaction.products[2]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.reactants[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[2]])
                    pairs.append([reaction.reactants[0], reaction.products[1]])
            elif reaction.products[2].containsLabeledAtom('*7'):
                if reaction.reactants[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[2]])
                elif reaction.reactants[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                    pairs.append([reaction.reactants[0], reaction.products[2]])
        elif reaction.isSurfaceReaction():
            # remove vacant active sites from consideration
            reactants = [sp for sp in reaction.reactants if not sp.isSurfaceSite()]
            products = [sp for sp in reaction.products if not sp.isSurfaceSite()]
            if len(reactants) == 1 or len(products) == 1:
                # When there is only one reactant (or one product), it is paired
                # with each of the products (reactants)
                for reactant in reactants:
                    for product in products:
                        pairs.append([reactant, product])
            elif self.label.lower() == 'surface_abstraction':
                # Hardcoding for surface abstraction: pair the reactant containing
                # *1 with the product containing *3 and vice versa
                assert len(reaction.reactants) == len(reaction.products) == 2
                if reaction.reactants[0].containsLabeledAtom('*1'):
                    if reaction.products[0].containsLabeledAtom('*3'):
                        pairs.append([reaction.reactants[0], reaction.products[0]])
                        pairs.append([reaction.reactants[1], reaction.products[1]])
                    elif reaction.products[1].containsLabeledAtom('*3'):
                        pairs.append([reaction.reactants[0], reaction.products[1]])
                        pairs.append([reaction.reactants[1], reaction.products[0]])
                elif reaction.reactants[1].containsLabeledAtom('*1'):
                    if reaction.products[1].containsLabeledAtom('*3'):
                        pairs.append([reaction.reactants[0], reaction.products[0]])
                        pairs.append([reaction.reactants[1], reaction.products[1]])
                    elif reaction.products[0].containsLabeledAtom('*3'):
                        pairs.append([reaction.reactants[0], reaction.products[1]])
                        pairs.append([reaction.reactants[1], reaction.products[0]])
        if not pairs:
            logging.debug('Preset mapping missing for determining reaction pairs for family {0!s}, '
                          'falling back to Reaction.generatePairs'.format(self.label))

        return pairs

    def getReactionTemplate(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, determine the most specific nodes in the tree that
        describe the reaction.
        """
        return self.groups.getReactionTemplate(reaction)

    def getKineticsForTemplate(self, template, degeneracy=1, method='rate rules'):
        """
        Return an estimate of the kinetics for a reaction with the given
        `template` and reaction-path `degeneracy`. There are two possible methods
        to use: 'group additivity' (new possible RMG-Py behavior) and 'rate rules' (old
        RMG-Java behavior, and default RMG-Py behavior).
        
        Returns a tuple (kinetics, entry):
        If it's estimated via 'rate rules' and an exact match is found in the tree,
        then the entry is returned as the second element of the tuple.
        But if an average is used, or the 'group additivity' method, then the tuple
        returned is (kinetics, None).
        """
        if method.lower() == 'group additivity':
            return self.estimateKineticsUsingGroupAdditivity(template, degeneracy), None
        elif method.lower() == 'rate rules':
            return self.estimateKineticsUsingRateRules(template, degeneracy)  # This returns kinetics and entry data
        else:
            raise ValueError('Invalid value "{0}" for method parameter; '
                             'should be "group additivity" or "rate rules".'.format(method))

    def getKineticsFromDepository(self, depository, reaction, template, degeneracy):
        """
        Search the given `depository` in this kinetics family for kinetics
        for the given `reaction`. Returns a list of all of the matching 
        kinetics, the corresponding entries, and ``True`` if the kinetics
        match the forward direction or ``False`` if they match the reverse
        direction.
        """
        kinetics_list = []
        entries = depository.entries.values()
        for entry in entries:
            if entry.item.isIsomorphic(reaction):
                kinetics_list.append(
                    [deepcopy(entry.data), entry, entry.item.isIsomorphic(reaction, eitherDirection=False)])
        for kinetics, entry, is_forward in kinetics_list:
            if kinetics is not None:
                kinetics.comment += "Matched reaction {0} {1} in {2}\nThis reaction matched rate rule {3}".format(
                    entry.index,
                    entry.label,
                    depository.label,
                    '[{0}]'.format(';'.join([g.label for g in template])))
                kinetics.comment += "\nfamily: {}".format(self.label)
        return kinetics_list

    def __selectBestKinetics(self, kineticsList):
        """
        For a given set of kinetics `kineticsList`, return the kinetics deemed
        to be the "best". This is determined to be the one with the lowest
        non-zero rank that occurs first (has the lowest index).
        """
        if any([x[1].rank == 0 for x in kineticsList]) and not all([x[1].rank == 0 for x in kineticsList]):
            kineticsList = [x for x in kineticsList if x[1].rank != 0]
        kineticsList.sort(key=lambda x: (x[1].rank, x[1].index))
        return kineticsList[0]

    def getKinetics(self, reaction, templateLabels, degeneracy=1, estimator='', returnAllKinetics=True):
        """
        Return the kinetics for the given `reaction` by searching the various
        depositories as well as generating a result using the user-specified `estimator`
        of either 'group additivity' or 'rate rules'.  Unlike
        the regular :meth:`getKinetics()` method, this returns a list of
        results, with each result comprising of

        1. the kinetics
        2. the source - this will be `None` if from a template estimate
        3. the entry  - this will be `None` if from a template estimate
        4. is_forward a boolean denoting whether the matched entry is in the same
           direction as the inputted reaction. This will always be True if using
           rates rules or group additivity. This can be `True` or `False` if using
           a depository

        If returnAllKinetics==False, only the first (best?) matching kinetics is returned.
        """
        kinetics_list = []

        depositories = self.depositories[:]

        template = self.retrieveTemplate(templateLabels)

        # Check the various depositories for kinetics
        for depository in depositories:
            kinetics_list0 = self.getKineticsFromDepository(depository, reaction, template, degeneracy)
            if len(kinetics_list0) > 0 and not returnAllKinetics:
                kinetics, entry, is_forward = self.__selectBestKinetics(kinetics_list0)
                return kinetics, depository, entry, is_forward
            else:
                for kinetics, entry, is_forward in kinetics_list0:
                    kinetics_list.append([kinetics, depository, entry, is_forward])

        # If estimator type of rate rules or group additivity is given, retrieve the kinetics. 
        if estimator:
            try:
                kinetics, entry = self.getKineticsForTemplate(template, degeneracy, method=estimator)
            except Exception:
                logging.error("Error getting kinetics for reaction {0!s}.\n{0!r}".format(reaction))
                raise

            if kinetics:
                if not returnAllKinetics:
                    return kinetics, estimator, entry, True
                kinetics_list.append([kinetics, estimator, entry, True])
        # If no estimation method was given, prioritize rate rule estimation. 
        # If returning all kinetics, add estimations from both rate rules and group additivity.
        else:
            try:
                kinetics, entry = self.getKineticsForTemplate(template, degeneracy, method='rate rules')
                if not returnAllKinetics:
                    return kinetics, 'rate rules', entry, True
                kinetics_list.append([kinetics, 'rate rules', entry, True])
            except KineticsError:
                # If kinetics were undeterminable for rate rules estimation, do nothing.
                pass

            try:
                kinetics2, entry2 = self.getKineticsForTemplate(template, degeneracy, method='group additivity')
                if not returnAllKinetics:
                    return kinetics2, 'group additivity', entry2, True
                kinetics_list.append([kinetics2, 'group additivity', entry2, True])
            except KineticsError:
                # If kinetics were undeterminable for group additivity estimation, do nothing.
                pass

        if not returnAllKinetics:
            raise UndeterminableKineticsError(reaction)

        return kinetics_list

    def estimateKineticsUsingGroupAdditivity(self, template, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using group additivity.
        
        Returns just the kinetics, or None.
        """
        warnings.warn("Group additivity is no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        # Start with the generic kinetics of the top-level nodes
        kinetics = None
        root = self.getRootTemplate()
        kinetics = self.getKineticsForTemplate(root)

        if kinetics is None:
            # raise UndeterminableKineticsError('Cannot determine group additivity kinetics estimate for '
            #                                   'template "{0}".'.format(','.join([e.label for e in template])))
            return None
        else:
            kinetics = kinetics[0]

        # Now add in more specific corrections if possible
        return self.groups.estimateKineticsUsingGroupAdditivity(template, kinetics, degeneracy)

    def estimateKineticsUsingRateRules(self, template, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using rate rules.
        
        Returns a tuple (kinetics, entry) where `entry` is the database
        entry used to determine the kinetics only if it is an exact match,
        and is None if some averaging or use of a parent node took place.
        """
        kinetics, entry = self.rules.estimateKinetics(template, degeneracy)

        return kinetics, entry

    def getReactionTemplateLabels(self, reaction):
        """
        Retrieve the template for the reaction and 
        return the corresponding labels for each of the 
        groups in the template.
        """
        template = self.getReactionTemplate(reaction)

        template_labels = []
        for entry in template:
            template_labels.append(entry.label)

        return template_labels

    def retrieveTemplate(self, templateLabels):
        """
        Reconstruct the groups associated with the 
        labels of the reaction template and 
        return a list.
        """
        template = []
        for label in templateLabels:
            template.append(self.groups.entries[label])

        return template

    def getLabeledReactantsAndProducts(self, reactants, products):
        """
        Given `reactants`, a list of :class:`Molecule` objects, and products, a list of 
        :class:`Molecule` objects, return two new lists of :class:`Molecule` objects with 
        atoms labeled: one for reactants, one for products. Returned molecules are totally 
        new entities in memory so input molecules `reactants` and `products` won't be affected.
        If RMG cannot find appropriate labels, (None, None) will be returned.
        """
        template = self.forwardTemplate
        reactants0 = [reactant.copy(deep=True) for reactant in reactants]

        if len(reactants0) == 1:
            molecule = reactants0[0]
            mappings = self.__matchReactantToTemplate(molecule, template.reactants[0].item)
            mappings = [[map0] for map0 in mappings]
            num_mappings = len(mappings)
            reactant_structures = [molecule]
        elif len(reactants0) == 2:
            molecule_a = reactants0[0]
            molecule_b = reactants0[1]
            # get mappings in forward direction
            mappings_a = self.__matchReactantToTemplate(molecule_a, template.reactants[0].item)
            mappings_b = self.__matchReactantToTemplate(molecule_b, template.reactants[1].item)
            mappings = list(itertools.product(mappings_a, mappings_b))
            # get mappings in the reverse direction
            mappings_a = self.__matchReactantToTemplate(molecule_a, template.reactants[1].item)
            mappings_b = self.__matchReactantToTemplate(molecule_b, template.reactants[0].item)
            mappings.extend(list(itertools.product(mappings_a, mappings_b)))

            reactant_structures = [molecule_a, molecule_b]
            num_mappings = len(mappings_a) * len(mappings_b)
        elif len(reactants0) == 3:
            molecule_a = reactants0[0]
            molecule_b = reactants0[1]
            molecule_c = reactants0[2]
            # Get mappings for all permutations of reactants
            mappings = []
            for order in itertools.permutations(range(3), 3):
                mappings_a = self.__matchReactantToTemplate(molecule_a, template.reactants[order[0]].item)
                mappings_b = self.__matchReactantToTemplate(molecule_b, template.reactants[order[1]].item)
                mappings_c = self.__matchReactantToTemplate(molecule_c, template.reactants[order[2]].item)
                mappings.extend(list(itertools.product(mappings_a, mappings_b, mappings_c)))

            reactant_structures = [molecule_a, molecule_b, molecule_c]
            num_mappings = len(mappings_a) * len(mappings_b) * len(mappings_c)
        else:
            raise IndexError('You have {0} reactants, which is unexpected!'.format(len(reactants)))

        for mapping in mappings:
            try:
                product_structures = self.__generateProductStructures(reactant_structures, mapping, forward=True)
            except ForbiddenStructureException:
                pass
            else:
                if product_structures is not None:
                    if same_species_lists(list(products), list(product_structures)):
                        return reactant_structures, product_structures
                    else:
                        continue

        # if there're some mapping available but cannot match the provided products
        # raise exception
        if num_mappings > 0:
            raise ActionError('Something wrong with products that RMG cannot find a match!')

        return None, None

    def addAtomLabelsForReaction(self, reaction, output_with_resonance=True):
        """
        Apply atom labels on a reaction using the appropriate atom labels from
        this reaction family.

        The reaction is modified in place containing species objects with the
        atoms labeled. If output_with_resonance is True, all resonance structures
        are generated with labels. If false, only the first resonance structure
        sucessfully able to map to the reaction is used. None is returned.
        """
        # make sure we start with reaction with species objects
        reaction.ensure_species(reactant_resonance=False, product_resonance=False)

        reactants = reaction.reactants
        products = reaction.products
        # ensure all species are independent references
        if len(reactants + products) > len(set([id(s) for s in reactants + products])):
            logging.debug('Copying reactants and products for reaction {} since they have '
                          'identical species references'.format(reaction))
            # not all species are independent
            reactants = [s.copy(deep=True) for s in reactants]
            products = [s.copy(deep=True) for s in products]

        # get all possible pairs of resonance structures
        reactant_pairs = list(itertools.product(*[s.molecule for s in reaction.reactants]))
        product_pairs = list(itertools.product(*[s.molecule for s in reaction.products]))

        labeled_reactants, labeled_products = None, None
        # go through each combination of possible pairs
        for reactant_pair, product_pair in itertools.product(reactant_pairs, product_pairs):
            try:
                # see if we obtain proper labeling
                labeled_reactants, labeled_products = self.getLabeledReactantsAndProducts(reactant_pair, product_pair)
                if labeled_reactants is not None:
                    break
            except ActionError:
                # must have gotten the wrong pair
                pass
        if labeled_reactants is None or labeled_products is None:
            raise ActionError("Could not find labeled reactants for reaction {} from "
                              "family {}.".format(reaction, self.label))

        # place the molecules in reaction's species object
        # this prevents overwriting of attributes of species objects by this method
        for index, species in enumerate(products):
            for labeled_molecule in labeled_products:
                if species.isIsomorphic(labeled_molecule):
                    species.molecule = [labeled_molecule]
                    reaction.products[index] = species
                    break
            else:
                raise ActionError('Could not find isomorphic molecule to fit the original product {} from '
                                  'reaction {}'.format(species, reaction))
        for index, species in enumerate(reactants):
            for labeled_molecule in labeled_reactants:
                if species.isIsomorphic(labeled_molecule):
                    species.molecule = [labeled_molecule]
                    reaction.reactants[index] = species
                    break
            else:
                raise ActionError('Could not find isomorphic molecule to fit the original reactant {} from '
                                  'reaction {}'.format(species, reaction))

        if output_with_resonance:
            # convert the molecules to species objects with resonance structures
            for species in reaction.reactants + reaction.products:
                species.generate_resonance_structures()

    def getw0(self, rxn):
        """
        calculates the w0 for Blower Masel kinetics by calculating wf (total bond energy of bonds formed)
        and wb (total bond energy of bonds broken) with w0 = (wf+wb)/2
        """
        mol = None
        aDict = {}
        for r in rxn.reactants:
            m = r.molecule[0]
            aDict.update(m.getLabeledAtoms())
            if mol:
                mol = mol.merge(m)
            else:
                mol = m.copy(deep=True)

        recipe = self.forwardRecipe.actions

        wb = 0.0
        wf = 0.0
        for act in recipe:

            if act[0] == 'BREAK_BOND':
                bd = mol.getBond(aDict[act[1]], aDict[act[3]])
                wb += bd.getBDE()
            elif act[0] == 'FORM_BOND':
                bd = Bond(aDict[act[1]], aDict[act[3]], act[2])
                wf += bd.getBDE()
            elif act[0] == 'CHANGE_BOND':
                bd1 = mol.getBond(aDict[act[1]], aDict[act[3]])

                if act[2] + bd1.order == 0.5:
                    mol2 = None
                    for r in rxn.products:
                        m = r.molecule[0]
                        if mol2:
                            mol2 = mol2.merge(m)
                        else:
                            mol2 = m.copy(deep=True)
                    bd2 = mol2.getBond(aDict[act[1]], aDict[act[3]])
                else:
                    bd2 = Bond(aDict[act[1]], aDict[act[3]], bd1.order + act[2])

                if bd2.order == 0:
                    bd2_bde = 0.0
                else:
                    bd2_bde = bd2.getBDE()
                bde_diff = bd2_bde - bd1.getBDE()
                if bde_diff > 0:
                    wf += abs(bde_diff)
                else:
                    wb += abs(bde_diff)

        return (wf + wb) / 2.0

    def getw0s(self, rxns):
        return list(map(self.getw0, rxns))

    def getTrainingDepository(self):
        """
        Returns the `training` depository from self.depositories
        """
        for depository in self.depositories:
            if depository.label.endswith('training'):
                return depository
        else:
            raise DatabaseError('Could not find training depository in family {0}.'.format(self.label))

    def addEntry(self, parent, grp, name):
        """
        Adds a group entry with parent parent
        group structure grp
        and group name name
        """
        ind = len(self.groups.entries) - 1
        entry = Entry(index=ind, label=name, item=grp, parent=parent)
        self.groups.entries[name] = entry
        self.rules.entries[name] = []
        if entry.parent:
            entry.parent.children.append(entry)

    def _splitReactions(self, rxns, newgrp):
        """
        divides the reactions in rxns between the new
        group structure newgrp and the old structure with 
        label oldlabel
        returns a list of reactions associated with the new group
        the list of reactions associated with the old group
        and a list of the indices of all of the reactions
        associated with the new group
        """
        new = []
        comp = []
        new_inds = []

        for i, rxn in enumerate(rxns):
            rmol = rxn.reactants[0].molecule[0]
            for r in rxn.reactants[1:]:
                rmol = rmol.merge(r.molecule[0])

            rmol.identifyRingMembership()

            if rmol.isSubgraphIsomorphic(newgrp, generateInitialMap=True, saveOrder=True):
                new.append(rxn)
                new_inds.append(i)
            else:
                comp.append(rxn)

        return new, comp, new_inds

    def evalExt(self, parent, ext, extname, templateRxnMap, obj=None, T=1000.0):
        """
        evaluates the objective function obj
        for the extension ext with name extname to the parent entry parent
        """
        rxns = templateRxnMap[parent.label]
        new, old, new_inds = self._splitReactions(rxns, ext)
        if len(new) == 0:
            return np.inf, False
        elif len(old) == 0:
            return np.inf, True
        else:
            if obj:
                ob, boo = getObjectiveFunction(new, old, obj, T=T)
            else:
                ob, boo = getObjectiveFunction(new, old, T=T)
            return ob, True

    def getExtensionEdge(self, parent, templateRxnMap, obj, T):
        """
        finds the set of all extension groups to parent such that
        1) the extension group divides the set of reactions under parent
        2) No generalization of the extension group divides the set of reactions under parent
                    
        We find this by generating all possible extensions of the initial group.  Extensions that split reactions are added
        to the list.  All extensions that do not split reactions and do not create bonds are ignored 
        (although those that match every reaction are labeled so we don't search them twice).  Those that match
        all reactions and involve bond creation undergo this process again.  
        
        Principle:  Say you have two elementary changes to a group ext1 and ext2 if applying ext1 and ext2 results in a 
        split at least one of ext1 and ext2 must result in a split
        
        Speed of this algorithm relies heavily on searching non bond creation dimensions once.
        """
        out_exts = [[]]
        grps = [parent.item]
        names = [parent.label]
        first_time = True

        n_splits = len(templateRxnMap[parent.label][0].reactants)

        while grps != []:
            grp = grps[-1]

            exts = grp.getExtensions(basename=names[-1], Nsplits=n_splits)

            reg_dict = dict()
            ext_inds = []
            for i, (grp2, grpc, name, typ, indc) in enumerate(exts):

                if typ != 'intNewBondExt' and typ != 'extNewBondExt' and (typ, indc) not in reg_dict.keys():
                    # first list is all extensions that match at least one reaction
                    # second is extensions that match all reactions
                    reg_dict[(typ, indc)] = ([], [])
                val, boo = self.evalExt(parent, grp2, name, templateRxnMap, obj, T)

                if val != np.inf:
                    out_exts[-1].append(exts[i])  # this extension splits reactions (optimization dim)
                    if typ == 'atomExt':
                        reg_dict[(typ, indc)][0].extend(grp2.atoms[indc[0]].atomType)
                    elif typ == 'elExt':
                        reg_dict[(typ, indc)][0].extend(grp2.atoms[indc[0]].radicalElectrons)
                    elif typ == 'bondExt':
                        reg_dict[(typ, indc)][0].extend(grp2.getBond(grp2.atoms[indc[0]], grp2.atoms[indc[1]]).order)

                elif boo:  # this extension matches all reactions (regularization dim)
                    if typ == 'intNewBondExt' or typ == 'extNewBondExt':
                        # these are bond formation extensions, we want to expand these until we get splits
                        ext_inds.append(i)
                    elif typ == 'atomExt':
                        reg_dict[(typ, indc)][0].extend(grp2.atoms[indc[0]].atomType)
                        reg_dict[(typ, indc)][1].extend(grp2.atoms[indc[0]].atomType)
                    elif typ == 'elExt':
                        reg_dict[(typ, indc)][0].extend(grp2.atoms[indc[0]].radicalElectrons)
                        reg_dict[(typ, indc)][1].extend(grp2.atoms[indc[0]].radicalElectrons)
                    elif typ == 'bondExt':
                        reg_dict[(typ, indc)][0].extend(grp2.getBond(grp2.atoms[indc[0]], grp2.atoms[indc[1]]).order)
                        reg_dict[(typ, indc)][1].extend(grp2.getBond(grp2.atoms[indc[0]], grp2.atoms[indc[1]]).order)
                    elif typ == 'ringExt':
                        reg_dict[(typ, indc)][1].append(True)
                else:
                    # this extension matches no reactions
                    if typ == 'ringExt':
                        reg_dict[(typ, indc)][0].append(False)
                        reg_dict[(typ, indc)][1].append(False)

            for typr, indcr in reg_dict.keys():  # have to label the regularization dimensions in all relevant groups
                reg_val = reg_dict[(typr, indcr)]

                if first_time and parent.children == []:
                    # parent
                    if typr != 'intNewBondExt' and typr != 'extNewBondExt':  # these dimensions should be regularized
                        if typr == 'atomExt':
                            grp.atoms[indcr[0]].reg_dim_atm = list(reg_val)
                        elif typr == 'elExt':
                            grp.atoms[indcr[0]].reg_dim_u = list(reg_val)
                        elif typr == 'ringExt':
                            grp.atoms[indcr[0]].reg_dim_r = list(reg_val)
                        elif typr == 'bondExt':
                            atms = grp.atoms
                            bd = grp.getBond(atms[indcr[0]], atms[indcr[1]])
                            bd.reg_dim = list(reg_val)

                # extensions being sent out
                if typr != 'intNewBondExt' and typr != 'extNewBondExt':  # these dimensions should be regularized
                    for grp2, grpc, name, typ, indc in out_exts[-1]:  # returned groups
                        if typr == 'atomExt':
                            grp2.atoms[indcr[0]].reg_dim_atm = list(reg_val)
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_atm = list(reg_val)
                        elif typr == 'elExt':
                            grp2.atoms[indcr[0]].reg_dim_u = list(reg_val)
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_u = list(reg_val)
                        elif typr == 'ringExt':
                            grp2.atoms[indcr[0]].reg_dim_r = list(reg_val)
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_r = list(reg_val)
                        elif typr == 'bondExt':
                            atms = grp2.atoms
                            bd = grp2.getBond(atms[indcr[0]], atms[indcr[1]])
                            bd.reg_dim = [list(set(bd.order) & set(reg_val[0])), list(set(bd.order) & set(reg_val[1]))]
                            if grpc:
                                atms = grpc.atoms
                                bd = grpc.getBond(atms[indcr[0]], atms[indcr[1]])
                                bd.reg_dim = [list(set(bd.order) & set(reg_val[0])),
                                              list(set(bd.order) & set(reg_val[1]))]

            # extensions being expanded
            for typr, indcr in reg_dict.keys():  # have to label the regularization dimensions in all relevant groups
                reg_val = reg_dict[(typr, indcr)]
                if typr != 'intNewBondExt' and typr != 'extNewBondExt':  # these dimensions should be regularized
                    for ind2 in ext_inds:  # groups for expansion
                        grp2, grpc, name, typ, indc = exts[ind2]
                        if typr == 'atomExt':
                            grp2.atoms[indcr[0]].reg_dim_atm = list(reg_val)
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_atm = list(reg_val)
                        elif typr == 'elExt':
                            grp2.atoms[indcr[0]].reg_dim_u = list(reg_val)
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_u = list(reg_val)
                        elif typr == 'ringExt':
                            grp2.atoms[indcr[0]].reg_dim_r = list(reg_val)
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_r = list(reg_val)
                        elif typr == 'bondExt':
                            atms = grp2.atoms
                            bd = grp2.getBond(atms[indcr[0]], atms[indcr[1]])
                            bd.reg_dim = [list(set(bd.order) & set(reg_val[0])), list(set(bd.order) & set(reg_val[1]))]
                            if grpc:
                                atms = grpc.atoms
                                bd = grpc.getBond(atms[indcr[0]], atms[indcr[1]])
                                bd.reg_dim = [list(set(bd.order) & set(reg_val[0])),
                                              list(set(bd.order) & set(reg_val[1]))]

            out_exts.append([])
            grps.pop()
            names.pop()

            for ind in ext_inds:  # collect the groups to be expanded
                grpr, grpcr, namer, typr, indcr = exts[ind]
                grps.append(grpr)
                names.append(namer)

            if first_time:
                first_time = False

        out = []
        # compile all of the valid extensions together
        # may be some duplicates here, but I don't think it's currently worth identifying them
        for x in out_exts:
            out.extend(x)

        return out

    def extendNode(self, parent, templateRxnMap, obj=None, T=1000.0, ):
        """
        Constructs an extension to the group parent based on evaluation 
        of the objective function obj
        """

        exts = self.getExtensionEdge(parent, templateRxnMap, obj=obj, T=T)

        if exts == []:  # should only occur when all reactions at this node are identical
            rs = templateRxnMap[parent.label]
            for q, rxn in enumerate(rs):
                for j in range(q):
                    if not same_species_lists(rxn.reactants, rs[j].reactants, generate_initial_map=True):
                        for p, atm in enumerate(parent.item.atoms):
                            if atm.reg_dim_atm[0] != atm.reg_dim_atm[1]:
                                logging.error('atom violation')
                                logging.error(atm.reg_dim_atm)
                                logging.error(parent.label)
                                logging.error('Regularization dimension suggest this node can be expanded, '
                                              'but extension generation has failed')
                            if atm.reg_dim_u[0] != atm.reg_dim_u[1]:
                                logging.error('radical violation')
                                logging.error(atm.reg_dim_u)
                                logging.error(parent.label)
                                logging.error('Regularization dimension suggest this node can be expanded, '
                                              'but extension generation has failed')
                        for p, bd in enumerate(parent.item.getAllEdges()):
                            if bd.reg_dim[0] != bd.reg_dim[1]:
                                logging.error('bond violation')
                                logging.error(bd.order)
                                logging.error(bd.reg_dim)
                                logging.error(parent.label)
                                logging.error('Regularization dimension suggest this node can be expanded, '
                                              'but extension generation has failed')

                        logging.error('split violation')
                        logging.error('parent')
                        logging.error(parent.item.toAdjacencyList())
                        for c, atm in enumerate(parent.item.atoms):
                            logging.error(c)
                            logging.error(atm.reg_dim_atm)
                            logging.error(atm.reg_dim_u)
                        logging.error("bonds:")
                        for bd in parent.item.getAllEdges():
                            ind1 = parent.item.atoms.index(bd.vertex1)
                            ind2 = parent.item.atoms.index(bd.vertex2)
                            logging.error(((ind1, ind2), bd.order, bd.reg_dim))
                        for rxn in rs:
                            for react in rxn.reactants:
                                logging.error(react.toAdjacencyList())
                        logging.error("Clearing Regularization Dimensions and Reattempting")  # this usually happens when node expansion breaks some symmetry
                        parent.item.clearRegDims()  # this almost always solves the problem
                        return True
            return False

        vals = []
        for grp, grpc, name, typ, einds in exts:
            val, boo = self.evalExt(parent, grp, name, templateRxnMap, obj, T)
            vals.append(val)

        min_val = min(vals)

        min_ind = vals.index(min_val)

        ext = exts[min_ind]

        extname = ext[2]

        if ext[3] == 'atomExt':
            ext[0].atoms[ext[4][0]].reg_dim_atm = [ext[0].atoms[ext[4][0]].atomType, ext[0].atoms[ext[4][0]].atomType]
        elif ext[3] == 'elExt':
            ext[0].atoms[ext[4][0]].reg_dim_u = [ext[0].atoms[ext[4][0]].radicalElectrons,
                                                 ext[0].atoms[ext[4][0]].radicalElectrons]

        self.addEntry(parent, ext[0], extname)

        complement = not ext[1] is None

        if complement:
            frags = extname.split('_')
            frags[-1] = 'N-' + frags[-1]
            cextname = ''
            for k in frags:
                cextname += k
                cextname += '_'
            cextname = cextname[:-1]

            self.addEntry(parent, ext[1], cextname)

        rxns = templateRxnMap[parent.label]

        new, left, new_inds = self._splitReactions(rxns, ext[0])

        comp_entries = []
        new_entries = []

        for i, entry in enumerate(templateRxnMap[parent.label]):
            if i in new_inds:
                new_entries.append(entry)
            else:
                comp_entries.append(entry)

        templateRxnMap[extname] = new_entries

        if complement:
            templateRxnMap[parent.label] = []
            templateRxnMap[cextname] = comp_entries
        else:
            templateRxnMap[parent.label] = comp_entries

        return True

    def generateTree(self, rxns=None, obj=None, thermoDatabase=None, T=1000.0, nprocs=1, minSplitableEntryNum=2,
                     minRxnsToSpawn=20, maxBatchSize=800, outlierFraction=0.02, stratumNum=8, maxRxnsToReoptNode=100):
        """
        Generate a tree by greedy optimization based on the objective function obj
        the optimization is done by iterating through every group and if the group has
        more than one training reaction associated with it a set of potential more specific extensions 
        are generated and the extension that optimizing the objective function combination is chosen 
        and the iteration starts over at the beginning
        
        additionally the tree structure is simplified on the fly by removing groups that have no kinetics data associated
        if their parent has no kinetics data associated and they either have only one child or
        have two children one of which has no kinetics data and no children
        (its parent becomes the parent of its only relevant child node)
        
        Args:
            rxns: List of reactions to generate tree from (if None pull the whole training set)
            obj: Object to expand tree from (if None uses top node)
            thermoDatabase: Thermodynamic database used for reversing training reactions 
            T: Temperature the tree is optimized for
            nprocs: Number of process for parallel tree generation 
            minSplitableEntryNum: the minimum number of splitable reactions at a node in order to spawn a new process solving that node
            minRxnsToSpawn: the minimum number of reactions at a node to spawn a new process solving that node
            maxBatchSize: the maximum number of reactions allowed in a batch, most batches will be this size the last will be smaller, 
                if the # of reactions < maxBatchSize the cascade algorithm is not used
            outlierFraction: Fraction of reactions that are fastest/slowest and will be automatically included in the first batch 
            stratumNum: Number of strata used in stratified sampling scheme 
            maxRxnsToReoptNode: Nodes with more matching reactions than this will not be pruned
        """
        if rxns is None:
            rxns = self.getTrainingSet(thermoDatabase=thermoDatabase, removeDegeneracy=True, estimateThermo=True,
                                       fixLabels=True, getReverse=True)

        if len(rxns) <= maxBatchSize:
            template_rxn_map = self.getReactionMatches(rxns=rxns, thermoDatabase=thermoDatabase, removeDegeneracy=True,
                                                       fixLabels=True, exactMatchesOnly=True, getReverse=True)
            self.makeTreeNodes(templateRxnMap=template_rxn_map, obj=obj, T=T, nprocs=nprocs - 1, depth=0,
                               minSplitableEntryNum=minSplitableEntryNum, minRxnsToSpawn=minRxnsToSpawn)
        else:
            random.seed(1)
            logging.error("dividing into batches")
            batches = self.getRxnBatches(rxns, T=T, maxBatchSize=maxBatchSize, outlierFraction=outlierFraction,
                                         stratumNum=stratumNum)
            logging.error([len(x) for x in batches])
            for i, batch in enumerate(batches):
                if i == 0:
                    rxns = batch
                else:
                    rxns += batch
                    logging.error("pruning tree")
                    self.pruneTree(rxns, thermoDatabase=thermoDatabase, maxRxnsToReoptNode=maxRxnsToReoptNode)
                logging.error("getting reaction matches")
                template_rxn_map = self.getReactionMatches(rxns=rxns, thermoDatabase=thermoDatabase, fixLabels=True,
                                                           exactMatchesOnly=True, getReverse=True)
                logging.error("building tree with {} rxns".format(len(rxns)))
                self.makeTreeNodes(templateRxnMap=template_rxn_map, obj=obj, T=T, nprocs=nprocs - 1, depth=0,
                                   minSplitableEntryNum=minSplitableEntryNum, minRxnsToSpawn=minRxnsToSpawn)

    def getRxnBatches(self, rxns, T=1000.0, maxBatchSize=800, outlierFraction=0.02, stratumNum=8):
        """
        Breaks reactions into batches based on a modified stratified sampling scheme
        Effectively:
        The top and bottom outlierFraction of all reactions are always included in the first batch
        The remaining reactions are ordered by the rate coefficients at T
        The list of reactions is then split into stratumNum similarly sized intervals
        batches sample equally from each interval, but randomly within each interval
        until they reach maxBatchSize reactions
        A list of lists of reactions containing the batches is returned
        """
        ks = np.array([rxn.kinetics.getRateCoefficient(T=T) for rxn in rxns])
        inds = np.argsort(ks)
        outlier_num = int(outlierFraction * len(ks) / 2)
        if outlier_num == 0:
            lowouts = []
            highouts = []
        else:
            lowouts = inds[:outlier_num].tolist()
            highouts = inds[-outlier_num:].tolist()
            inds = inds[outlier_num:-outlier_num]
        interval_length = int(len(inds) / stratumNum)
        strata = []
        for i in range(stratumNum):
            if i == 0:
                temp = inds[:interval_length].tolist()
                random.shuffle(temp)
                strata.append(temp)
            elif i == stratumNum - 1:
                temp = inds[interval_length * i:].tolist()
                random.shuffle(temp)
                strata.append(temp)
            else:
                temp = inds[interval_length * i:interval_length * (i + 1)].tolist()
                random.shuffle(temp)
                strata.append(temp)

        first_batch_strata_num = maxBatchSize - outlier_num
        batches = [highouts + lowouts]
        bind = 0
        while any([len(stratum) != 0 for stratum in strata]):
            for stratum in strata:
                if stratum != []:
                    batches[bind].append(stratum.pop())
                    if len(batches[bind]) >= maxBatchSize:
                        bind += 1
                        batches.append([])

        rxns = np.array(rxns)
        batches = [rxns[inds].tolist() for inds in batches if len(inds) > 0]

        return batches

    def pruneTree(self, rxns, thermoDatabase=None, maxRxnsToReoptNode=100, fixLabels=True, exactMatchesOnly=True,
                  getReverse=True):
        """
        Remove nodes that have less than maxRxnToReoptNode reactions that match
        and clear the regularization dimensions of their parent
        This is used to remove smaller easier to optimize and more likely to change nodes
        before adding a new batch in cascade model generation
        """
        template_rxn_map = self.getReactionMatches(rxns=rxns, thermoDatabase=thermoDatabase, fixLabels=fixLabels,
                                                   exactMatchesOnly=False, getReverse=getReverse)
        for key, item in template_rxn_map.items():
            entry = self.groups.entries[key]
            parent = entry.parent
            if parent and len(template_rxn_map[parent.label]) < maxRxnsToReoptNode:
                parent.children.remove(entry)
                del self.groups.entries[key]
                parent.item.clearRegDims()

    def makeTreeNodes(self, templateRxnMap=None, obj=None, T=1000.0, nprocs=0, depth=0, minSplitableEntryNum=2,
                      minRxnsToSpawn=20):

        if depth > 0:
            root = self.groups.entries[list(templateRxnMap.keys())[0]]
        else:
            for entry in self.groups.entries.values():  # find the root entry for this branch
                if entry.index != -1:
                    root = entry
                    break
            while root.parent is not None:
                root = root.parent

        psize = float(len(templateRxnMap[root.label]))

        mult_completed_nodes = []  # nodes containing multiple identical training reactions
        boo = True  # if the for loop doesn't break becomes false and the while loop terminates
        active_procs = []
        active_conns = []
        active_proc_num = []
        proc_names = []
        free_procs = nprocs
        extra_entries = []

        while boo:
            remove_inds = []
            for k, p in enumerate(active_procs):  # check if any processes have finished
                if active_conns[k].poll():
                    new_entries = self._absorbProcess(p, active_conns[k], proc_names[k])
                    extra_entries += new_entries
                    remove_inds.append(k)

            remove_inds.reverse()
            for ind in remove_inds:  # remove finished process objects
                free_procs += active_proc_num[ind]
                del active_proc_num[ind]
                del active_procs[ind]
                del active_conns[ind]
                del proc_names[ind]

            splitable_entry_num = 0
            for label, items in templateRxnMap.items():  # figure out how many splitable objects there are
                entry = self.groups.entries[label]
                if len(items) > 1 and entry not in mult_completed_nodes:
                    splitable_entry_num += 1

            for label in templateRxnMap.keys():
                entry = self.groups.entries[label]
                if not isinstance(entry.item, Group):  # skip logic nodes
                    continue
                if psize == 0.0:
                    continue
                if entry.index != -1 and len(templateRxnMap[entry.label]) > 1 and entry not in mult_completed_nodes:
                    if (free_procs > 0 and splitable_entry_num > minSplitableEntryNum and
                            len(templateRxnMap[entry.label]) > minRxnsToSpawn):
                        procs_out = int(len(templateRxnMap[entry.label]) / psize * free_procs)
                        free_procs -= procs_out
                        assert free_procs >= 0
                        conn, p, name = _spawnTreeProcess(family=self,
                                                          templateRxnMap={entry.label: templateRxnMap[entry.label]},
                                                          obj=obj, T=T, nprocs=procs_out - 1, depth=depth,
                                                          minSplitableEntryNum=minSplitableEntryNum,
                                                          minRxnsToSpawn=minRxnsToSpawn)
                        active_procs.append(p)
                        active_conns.append(conn)
                        proc_names.append(name)
                        active_proc_num.append(procs_out)
                        L = entry.label
                        self.groups.entries[L].parent.children.remove(entry)
                        del templateRxnMap[L]  # prevents this process from recreating work done by another process
                        del self.groups.entries[L]

                        splitable_entry_num -= 1
                        continue
                    boo2 = self.extendNode(entry, templateRxnMap, obj, T)
                    if boo2:  # extended node so restart while loop
                        break
                    else:  # no extensions could be generated since all reactions were identical
                        mult_completed_nodes.append(entry)
            else:
                if len(active_procs) == 0:
                    boo = False

            # fix indicies
            iters = 0
            for entry in self.groups.entries.values():
                if entry.index != -1:
                    entry.index = iters
                    iters += 1

        # add the entries generated on other processors
        index = max([ent.index for ent in self.groups.entries.values()]) + 1
        for item in extra_entries:
            if item.label in self.groups.entries.keys():
                continue
            item.index = index
            index += 1
            self.groups.entries[item.label] = item

        for label, entry in self.groups.entries.items():
            if entry.index != -1 and entry.parent is None and entry.label != root.label:
                pname = "_".join(label.split('_')[:-1])
                while pname not in self.groups.entries.keys():
                    pname = "_".join(label.split('_')[:-1])
                entry.parent = self.groups.entries[pname]
                entry.parent.children.append(entry)

        return

    def _absorbProcess(self, p, conn, name):
        try:
            grps = conn.recv()
            p.terminate()
        except Exception as e:
            logging.error('failed to absorb process {}'.format(name))
            raise e
        return grps

    def makeBMRulesFromTemplateRxnMap(self, templateRxnMap, nprocs=1, Tref=1000.0, fmax=1.0e5):

        rule_keys = self.rules.entries.keys()
        for entry in self.groups.entries.values():
            if entry.label not in rule_keys:
                self.rules.entries[entry.label] = []

        index = max([e.index for e in self.rules.getEntries()] or [0]) + 1

        entries = self.groups.entries.values()
        rxnlists = [(templateRxnMap[entry.label], entry.label)
                    if entry.label in templateRxnMap.keys() else [] for entry in entries]
        inputs = np.array([(self.forwardRecipe.actions, rxns, Tref, fmax, label, [r.rank for r in rxns])
                           for rxns, label in rxnlists])

        inds = np.arange(len(inputs))
        np.random.shuffle(inds)  # want to parallelize in random order
        inds = inds.tolist()
        revinds = [inds.index(x) for x in np.arange(len(inputs))]

        pool = mp.Pool(nprocs)

        kinetics_list = np.array(pool.map(_makeRule, inputs[inds]))
        kinetics_list = kinetics_list[revinds]  # fix order

        for i, kinetics in enumerate(kinetics_list):
            if kinetics is not None:
                entry = entries[i]
                std = kinetics.uncertainty.getExpectedLogUncertainty() / 0.398  # expected uncertainty is std * 0.398
                st = "BM rule fitted to {0} training reactions at node {1}".format(len(rxnlists[i]), entry.label)
                st += "\nTotal Standard Deviation in ln(k): {0}".format(std)
                new_entry = Entry(
                    index=index,
                    label=entry.label,
                    item=self.forwardTemplate,
                    data=kinetics,
                    rank=11,
                    reference=None,
                    shortDesc=st,
                    longDesc=st,
                )
                new_entry.data.comment = st

                self.rules.entries[entry.label].append(new_entry)

                index += 1

    def crossValidate(self, folds=5, templateRxnMap=None, testRxnInds=None, T=1000.0, iters=0, random_state=1):
        """
        Perform K-fold cross validation on an automatically generated tree at temperature T
        after finding an appropriate node for kinetics estimation it will move up the tree
        iters times.  
        Returns a dictionary mapping {rxn:Ln(k_Est/k_Train)}
        """

        if templateRxnMap is None:
            templateRxnMap = self.getReactionMatches(removeDegeneracy=True, getReverse=True, fixLabels=True)

        rxns = np.array(templateRxnMap['Root'])

        if testRxnInds is None:
            if folds == 0:
                folds = len(rxns)

            kf = KFold(folds, shuffle=True, random_state=random_state)
            kfsplits = kf.split(rxns)
        else:
            kfsplits = [([0, ], [0, ])]

        errors = {}
        uncertainties = {}

        for train_index, test_index in kfsplits:

            if testRxnInds is None:
                rxns_test = rxns[test_index]
            else:
                rxns_test = rxns[testRxnInds]

            for rxn in rxns_test:

                krxn = rxn.kinetics.getRateCoefficient(T)

                entry = self.getRootTemplate()[0]

                boo = True
                while boo:  # find the entry it matches
                    for child in entry.children:
                        rs = templateRxnMap[child.label]
                        if rxn in rs:
                            entry = child
                            break
                    else:
                        boo = False

                while entry.parent and len(set(templateRxnMap[entry.label]) - set(rxns_test)) <= 1:
                    if entry.parent:
                        entry = entry.parent

                for q in range(iters):
                    if entry.parent:
                        entry = entry.parent

                uncertainties[rxn] = self.rules.entries[entry.label][0].data.uncertainty

                L = list(set(templateRxnMap[entry.label]) - set(rxns_test))

                if L != []:
                    kinetics = ArrheniusBM().fitToReactions(L, recipe=self.forwardRecipe.actions)
                    kinetics = kinetics.toArrhenius(rxn.getEnthalpyOfReaction(T))
                    k = kinetics.getRateCoefficient(T)
                    errors[rxn] = np.log(k / krxn)
                else:
                    raise ValueError('only one piece of kinetics information in the tree?')

        return errors, uncertainties

    def crossValidateOld(self, folds=5, T=1000.0, random_state=1, estimator='rate rules', thermoDatabase=None):
        """
        Perform K-fold cross validation on an automatically generated tree at temperature T
        Returns a dictionary mapping {rxn:Ln(k_Est/k_Train)}
        """
        errors = {}
        rxns = np.array(self.getTrainingSet(removeDegeneracy=True))

        if folds == 0:
            folds = len(rxns)

        kf = KFold(folds, shuffle=True, random_state=random_state)

        if thermoDatabase is None:
            from rmgpy.data.rmg import getDB
            tdb = getDB('thermo')
        else:
            tdb = thermoDatabase

        for train_index, test_index in kf.split(rxns):

            self.rules.entries = {}  # clear rules each iteration

            self.addKineticsRulesFromTrainingSet(trainIndices=train_index, thermoDatabase=tdb)
            self.fillKineticsRulesByAveragingUp()
            rxns_test = rxns[test_index]

            for rxn in rxns_test:

                krxn = rxn.kinetics.getRateCoefficient(T)

                template_labels = self.getReactionTemplateLabels(rxn)
                template = self.retrieveTemplate(template_labels)
                if estimator == 'rate rules':
                    kinetics, entry = self.estimateKineticsUsingRateRules(template, degeneracy=1)
                elif estimator == 'group additivity':
                    kinetics = self.estimateKineticsUsingGroupAdditivity(template, degeneracy=1)
                else:
                    raise ValueError('{0} is not a valid value for input `estimator`'.format(estimator))

                k = kinetics.getRateCoefficient(T)

                errors[rxn] = np.log(k / krxn)

        return errors

    def simpleRegularization(self, node, templateRxnMap, test=True):
        """
        Simplest regularization algorithm
        All nodes are made as specific as their descendant reactions
        Training reactions are assumed to not generalize 
        For example if an particular atom at a node is Oxygen for all of its
        descendent reactions a reaction where it is Sulfur will never hit that node
        unless it is the top node even if the tree did not split on the identity 
        of that atom
        
        The test option to this function determines whether or not the reactions 
        under a node match the extended group before adding an extension. 
        If the test fails the extension is skipped. 
        
        In general test=True is needed if the cascade algorithm was used 
        to generate the tree and test=False is ok if the cascade algorithm
        wasn't used. 
        """

        for child in node.children:
            self.simpleRegularization(child, templateRxnMap)

        grp = node.item
        rxns = templateRxnMap[node.label]

        R = ['H', 'C', 'N', 'O', 'Si', 'S', 'Cl']  # set of possible R elements/atoms
        R = [atomTypes[x] for x in R]

        RnH = R[:]
        RnH.remove(atomTypes['H'])

        Run = [0, 1, 2, 3]

        atm_dict = {'R': R, 'R!H': RnH}

        if isinstance(node.item, Group):
            indistinguishable = []
            for i, atm1 in enumerate(grp.atoms):

                skip = False
                if node.children == []:  # if the atoms or bonds are graphically indistinguishable don't regularize
                    bdpairs = {(atm, tuple(bd.order)) for atm, bd in atm1.bonds.items()}
                    for atm2 in grp.atoms:
                        if atm1 is not atm2 and atm1.atomType == atm2.atomType and len(atm1.bonds) == len(atm2.bonds):
                            bdpairs2 = {(atm, tuple(bd.order)) for atm, bd in atm2.bonds.items()}
                            if bdpairs == bdpairs2:
                                skip = True
                                indistinguishable.append(i)

                if not skip and atm1.reg_dim_atm[1] != [] and set(atm1.reg_dim_atm[1]) != set(atm1.atomType):
                    atyp = atm1.atomType
                    if len(atyp) == 1 and atyp[0] in R:
                        pass
                    else:
                        if len(atyp) == 1 and atyp[0].label in atm_dict.keys():
                            atyp = atm_dict[atyp[0].label]

                        vals = list(set(atyp) & set(atm1.reg_dim_atm[1]))
                        assert vals != [], 'cannot regularize to empty'
                        if all([set(child.item.atoms[i].atomType) <= set(vals) for child in node.children]):
                            if not test:
                                atm1.atomType = vals
                            else:
                                oldvals = atm1.atomType
                                atm1.atomType = vals
                                if not self.rxnsMatchNode(node, rxns):
                                    atm1.atomType = oldvals

                if not skip and atm1.reg_dim_u[1] != [] and set(atm1.reg_dim_u[1]) != set(atm1.radicalElectrons):
                    if len(atm1.radicalElectrons) == 1:
                        pass
                    else:
                        relist = atm1.radicalElectrons
                        if relist == []:
                            relist = Run
                        vals = list(set(relist) & set(atm1.reg_dim_u[1]))
                        assert vals != [], 'cannot regularize to empty'

                        if all([set(child.item.atoms[i].radicalElectrons) <= set(vals)
                                if child.item.atoms[i].radicalElectrons != [] else False for child in node.children]):
                            if not test:
                                atm1.radicalElectrons = vals
                            else:
                                oldvals = atm1.radicalElectrons
                                atm1.radicalElectrons = vals
                                if not self.rxnsMatchNode(node, rxns):
                                    atm1.radicalElectrons = oldvals

                if (not skip and atm1.reg_dim_r[1] != [] and
                        ('inRing' not in atm1.props.keys() or atm1.reg_dim_r[1][0] != atm1.props['inRing'])):
                    if 'inRing' not in atm1.props.keys():
                        if (all(['inRing' in child.item.atoms[i].props.keys() for child in node.children]) and
                                all([child.item.atoms[i].props['inRing'] == atm1.reg_dim_r[1] for child in node.children])):
                            if not test:
                                atm1.props['inRing'] = atm1.reg_dim_r[1][0]
                            else:
                                if 'inRing' in atm1.props.keys():
                                    oldvals = atm1.props['inRing']
                                else:
                                    oldvals = None
                                atm1.props['inRing'] = atm1.reg_dim_r[1][0]
                                if not self.rxnsMatchNode(node, rxns):
                                    if oldvals:
                                        atm1.props['inRing'] = oldvals
                                    else:
                                        del atm1.props['inRing']
                if not skip:
                    for j, atm2 in enumerate(grp.atoms[:i]):
                        if j in indistinguishable:  # skip graphically indistinguishable atoms
                            continue
                        if grp.hasBond(atm1, atm2):
                            bd = grp.getBond(atm1, atm2)
                            if len(bd.order) == 1:
                                pass
                            else:
                                vals = list(set(bd.order) & set(bd.reg_dim[1]))
                                if vals != [] and all([set(child.item.getBond(child.item.atoms[i], child.item.atoms[j]).order) <= set(vals) for child in node.children]):
                                    if not test:
                                        bd.order = vals
                                    else:
                                        oldvals = bd.order
                                        bd.order = vals
                                        if not self.rxnsMatchNode(node, rxns):
                                            bd.order = oldvals

    def regularize(self, regularization=simpleRegularization, keepRoot=True, thermoDatabase=None, templateRxnMap=None,
                   rxns=None):
        """
        Regularizes the tree according to the regularization function regularization
        """
        if templateRxnMap is None:
            if rxns is None:
                templateRxnMap = self.getReactionMatches(thermoDatabase=thermoDatabase, removeDegeneracy=True,
                                                         getReverse=True, exactMatchesOnly=False, fixLabels=True)
            else:
                templateRxnMap = self.getReactionMatches(rxns=rxns, thermoDatabase=thermoDatabase,
                                                         removeDegeneracy=True, getReverse=True, exactMatchesOnly=False,
                                                         fixLabels=True)

        if keepRoot:
            for child in self.getRootTemplate()[0].children:  # don't regularize the root
                regularization(self, child, templateRxnMap)
        else:
            regularization(self, self.getRootTemplate()[0], templateRxnMap)

    def checkTree(self, entry=None):
        if entry is None:
            entry = self.getRootTemplate()[0]
        for child in entry.children:
            if not child.item.isSubgraphIsomorphic(entry.item, generateInitialMap=True, saveOrder=True):
                logging.error('child: ')
                logging.error(child.label)
                logging.error(child.item.toAdjacencyList())
                logging.error('parent: ')
                logging.error(entry.label)
                logging.error(entry.item.toAdjacencyList())
                raise ValueError('Child not subgraph isomorphic to parent')
            self.checkTree(child)
        for entry in self.groups.entries.values():
            if entry.index == -1:
                continue
            parent = entry
            while parent.parent is not None:
                parent = parent.parent
            assert parent.label == 'Root', parent.label

    def makeTree(self, obj=None, regularization=simpleRegularization, thermoDatabase=None, T=1000.0):
        """
        generates tree structure and then generates rules for the tree
        """
        self.generateTree(obj=obj, thermoDatabase=thermoDatabase, T=T)
        self.regularize(regularization=regularization)
        template_rxn_map = self.getReactionMatches(thermoDatabase=thermoDatabase,
                                                   removeDegeneracy=True, getReverse=True)
        self.makeBMRulesFromTemplateRxnMap(template_rxn_map)
        self.checkTree()

    def cleanTreeRules(self):
        self.rules.entries = OrderedDict()
        self.rules.entries['Root'] = []

    def cleanTreeGroups(self, thermoDatabase=None):
        """
        clears groups and rules in the tree, generates an appropriate
        root group to start from and then reads training reactions
        Note this only works if a single top node (not a logic node)
        can be generated
        """
        # find the starting node
        grp = None

        rtmps = self.getRootTemplate()

        if not isinstance(rtmps[0].item, Group):
            raise ValueError('each tree top node must be a group not a logic node to prepare the tree automatically')

        for ent in rtmps:
            if grp is None:
                grp = ent.item
            else:
                if any([isinstance(x, list) for x in ent.item.getLabeledAtoms().values()]):
                    grp = grp.mergeGroups(ent.item, keepIdenticalLabels=True)
                else:
                    grp = grp.mergeGroups(ent.item)

        # clear everything
        self.groups.entries = {x.label: x for x in self.groups.entries.values() if x.index == -1}

        # add the starting node
        self.addEntry(None, grp, 'Root')
        self.groups.entries['Root'].index = 1
        self.groups.top = [self.groups.entries['Root']]
        self.forwardTemplate.reactants = [self.groups.entries['Root']]

        return

    def cleanTree(self, thermoDatabase=None):
        self.cleanTreeRules()
        self.cleanTreeGroups(thermoDatabase=thermoDatabase)

    def saveGeneratedTree(self, path=None):
        """
        clears the rules and saves the family to its 
        current location in database
        """
        if path is None:
            path = settings['database.directory']
            path = os.path.join(path, 'kinetics', 'families', self.label)

        self.save(path)

    def getTrainingSet(self, thermoDatabase=None, removeDegeneracy=False, estimateThermo=True, fixLabels=False,
                       getReverse=False):
        """
        retrieves all reactions in the training set, assigns thermo to the species objects
        reverses reactions as necessary so that all reactions are in the forward direction
        and returns the resulting list of reactions in the forward direction with thermo 
        assigned
        """

        def getLabelFixedMol(mol, root_labels):
            nmol = mol.copy(deep=True)
            for atm in nmol.atoms:
                if atm.label not in root_labels:
                    atm.label = ''
            return nmol

        def fixLabelsMol(mol, root_labels):
            for atm in mol.atoms:
                if atm.label not in root_labels:
                    atm.label = ''

        if self.ownReverse and getReverse:
            rev_rxns = []
            rkeys = list(self.reverseMap.keys())
            reverse_map = self.reverseMap

        if estimateThermo:
            if thermoDatabase is None:
                from rmgpy.data.rmg import getDB
                tdb = getDB('thermo')
            else:
                tdb = thermoDatabase

        try:
            dep = self.getTrainingDepository()
        except:
            logging.info('Could not find training depository in family {0}.'.format(self.label))
            logging.info('Must be because you turned off the training depository.')
            return

        rxns = deepcopy([i.item for i in dep.entries.values()])
        entries = deepcopy([i for i in dep.entries.values()])

        roots = [x.item for x in self.getRootTemplate()]
        root = None
        for r in roots:
            if root:
                root = root.mergeGroups(r)
            else:
                root = deepcopy(r)

        root_labels = [x.label for x in root.atoms if x.label != '']

        for i, r in enumerate(entries):
            if estimateThermo:
                for j, react in enumerate(r.item.reactants):
                    if rxns[i].reactants[j].thermo is None:
                        rxns[i].reactants[j].thermo = tdb.getThermoData(react)

                for j, react in enumerate(r.item.products):
                    if rxns[i].products[j].thermo is None:
                        rxns[i].products[j].thermo = tdb.getThermoData(react)

            rxns[i].kinetics = r.data
            rxns[i].rank = r.rank

            if removeDegeneracy:  # adjust for degeneracy
                rxns[i].kinetics.A.value_si /= rxns[i].degeneracy

            mol = None
            for react in rxns[i].reactants:
                if fixLabels:
                    fixLabelsMol(react.molecule[0], root_labels)
                if mol:
                    mol = mol.merge(react.molecule[0])
                else:
                    mol = deepcopy(react.molecule[0])

            if fixLabels:
                for prod in rxns[i].products:
                    fixLabelsMol(prod.molecule[0], root_labels)
                for atm in mol.atoms:
                    if atm.label not in root_labels:
                        atm.label = ''

            if (mol.isSubgraphIsomorphic(root, generateInitialMap=True) or
                    (not fixLabels and
                     getLabelFixedMol(mol, root_labels).isSubgraphIsomorphic(root, generateInitialMap=True))):
                rxns[i].is_forward = True
                if self.ownReverse and getReverse:
                    mol = None
                    for react in rxns[i].products:
                        if mol:
                            mol = mol.merge(react.molecule[0])
                        else:
                            mol = deepcopy(react.molecule[0])

                    if (mol.isSubgraphIsomorphic(root, generateInitialMap=True) or
                            (not fixLabels and
                             getLabelFixedMol(mol, root_labels).isSubgraphIsomorphic(root, generateInitialMap=True))):
                        # try product structures
                        products = [Species(molecule=[getLabelFixedMol(x.molecule[0], root_labels)], thermo=x.thermo)
                                    for x in rxns[i].products]
                    else:
                        products = self.applyRecipe([s.molecule[0] for s in rxns[i].reactants], forward=True)
                        products = [Species(molecule=[p]) for p in products]

                    prodmol = None
                    for react in rxns[i].products:
                        if prodmol:
                            prodmol = prodmol.merge(react.molecule[0])
                        else:
                            prodmol = deepcopy(react.molecule[0])

                    if not prodmol.isSubgraphIsomorphic(root, generateInitialMap=True):
                        mol = None
                        for react in products:
                            if mol:
                                mol = mol.merge(react.molecule[0])
                            else:
                                mol = deepcopy(react.molecule[0])
                        if not mol.isSubgraphIsomorphic(root, generateInitialMap=True):
                            for p in products:
                                for atm in p.molecule[0].atoms:
                                    if atm.label in rkeys:
                                        atm.label = reverse_map[atm.label]

                    reacts = [Species(molecule=[getLabelFixedMol(x.molecule[0], root_labels)], thermo=x.thermo)
                              for x in rxns[i].reactants]
                    rrev = Reaction(reactants=products, products=reacts,
                                    kinetics=rxns[i].generateReverseRateCoefficient(), rank=rxns[i].rank)
                    rrev.is_forward = False

                    if estimateThermo:
                        for r in rrev.reactants:
                            if r.thermo is None:
                                r.thermo = tdb.getThermoData(deepcopy(r))

                    rev_rxns.append(rrev)

                continue
            else:
                if self.ownReverse:
                    logging.error("rxn")
                    logging.error(str(rxns[i]))
                    logging.error("root")
                    logging.error(root.toAdjacencyList())
                    logging.error("mol")
                    logging.error(mol.toAdjacencyList())
                    raise ValueError("couldn't match reaction")

                mol = None
                for react in rxns[i].products:
                    if mol:
                        mol = mol.merge(react.molecule[0])
                    else:
                        mol = deepcopy(react.molecule[0])

                if (mol.isSubgraphIsomorphic(root, generateInitialMap=True) or
                        (not fixLabels and
                         getLabelFixedMol(mol, root_labels).isSubgraphIsomorphic(root, generateInitialMap=True))):  # try product structures
                    products = [Species(molecule=[getLabelFixedMol(x.molecule[0], root_labels)], thermo=x.thermo)
                                for x in rxns[i].products]
                else:
                    products = self.applyRecipe([s.molecule[0] for s in rxns[i].reactants], forward=True)
                    products = [Species(molecule=[p]) for p in products]

                rrev = Reaction(reactants=products, products=rxns[i].reactants,
                                kinetics=rxns[i].generateReverseRateCoefficient(), rank=rxns[i].rank)

                rrev.is_forward = False

                if estimateThermo:
                    for r in rrev.reactants:
                        if r.thermo is None:
                            r.thermo = tdb.getThermoData(deepcopy(r))
                rxns[i] = rrev

        if self.ownReverse and getReverse:
            return rxns + rev_rxns
        else:
            return rxns

    def getReactionMatches(self, rxns=None, thermoDatabase=None, removeDegeneracy=False, estimateThermo=True,
                           fixLabels=False, exactMatchesOnly=False, getReverse=False):
        """
        returns a dictionary mapping for each entry in the tree:  
        (entry.label,entry.item) : list of all training reactions (or the list given) that match that entry
        """
        if rxns is None:
            rxns = self.getTrainingSet(thermoDatabase=thermoDatabase, removeDegeneracy=removeDegeneracy,
                                       estimateThermo=estimateThermo, fixLabels=fixLabels, getReverse=getReverse)

        entries = self.groups.entries

        assert len(set(entries.keys())) == len(entries.keys()), 'there are duplicate indices in family.group.entries'

        rxn_lists = {entry.label: [] for entry in entries.values()}

        root = self.getRootTemplate()[0]

        for rxn in rxns:
            mol = None
            for r in rxn.reactants:
                if mol is None:
                    mol = deepcopy(r.molecule[0])
                else:
                    mol = mol.merge(r.molecule[0])
            try:
                flag = not self.isEntryMatch(mol, root, resonance=True)
            except:
                flag = not self.isEntryMatch(mol, root, resonance=False)

            if flag:
                logging.error(root.item.toAdjacencyList())
                logging.error(mol.toAdjacencyList())
                for r in rxn.reactants:
                    logging.error(r.molecule[0].toAdjacencyList())
                for r in rxn.products:
                    logging.error(r.molecule[0].toAdjacencyList())
                raise ValueError('reaction: {0} does not match root template in family {1}'.format(rxn, self.label))

            rxn_lists[root.label].append(rxn)

            entry = root

            while entry.children != []:
                for child in entry.children:
                    if self.isEntryMatch(mol, child, resonance=False):
                        entry = child
                        rxn_lists[child.label].append(rxn)
                        break
                else:
                    break

        if exactMatchesOnly:
            new_lists = dict()
            for key, rs in rxn_lists.items():
                newrs = set(rs)
                for child in self.groups.entries[key].children:
                    newrs -= set(rxn_lists[child.label])
                new_lists[key] = list(newrs)
            rxn_lists = new_lists

        return rxn_lists

    def isEntryMatch(self, mol, entry, resonance=True):
        """
        determines if the labeled molecule object of reactants matches the entry entry
        """
        if isinstance(entry.item, Group):
            if resonance:
                structs = mol.generate_resonance_structures()
            else:
                structs = [mol]
            return any([mol.isSubgraphIsomorphic(entry.item, generateInitialMap=True) for mol in structs])
        elif isinstance(entry.item, LogicOr):
            return any([self.isEntryMatch(mol, self.groups.entries[c], resonance=resonance)
                        for c in entry.item.components])

    def rxnsMatchNode(self, node, rxns):
        for rxn in rxns:
            mol = None
            for r in rxn.reactants:
                if mol is None:
                    mol = deepcopy(r.molecule[0])
                else:
                    mol = mol.merge(r.molecule[0])

            if not self.isEntryMatch(mol, node, resonance=False):
                return False

        return True

    def retrieveOriginalEntry(self, templateLabel):
        """
        Retrieves the original entry, be it a rule or training reaction, given
        the template label in the form 'group1;group2' or 'group1;group2;group3'
        
        Returns tuple in the form
        (RateRuleEntry, TrainingReactionEntry)
        
        Where the TrainingReactionEntry is only present if it comes from a training reaction
        """
        template_labels = templateLabel.split()[-1].split(';')
        template = self.retrieveTemplate(template_labels)
        rule = self.getRateRule(template)
        if 'From training reaction' in rule.data.comment:
            training_index = int(rule.data.comment.split()[3])
            training_depository = self.getTrainingDepository()
            return rule, training_depository.entries[training_index]
        else:
            return rule, None

    def getSourcesForTemplate(self, template):
        """
        Returns the set of rate rules and training reactions used to average this `template`.  Note that the tree must be
        averaged with verbose=True for this to work.
        
        Returns a tuple of
        rules, training
        
        where rules are a list of tuples containing 
        the [(original_entry, weight_used_in_average), ... ]
        
        and training is a list of tuples containing
        the [(rate_rule_entry, training_reaction_entry, weight_used_in_average),...]
        """

        def assignWeightsToEntries(entry_nested_list, weighted_entries, n=1):
            """
            Assign weights to an average of average nested list. Where n is the
            number of values being averaged recursively.  
            """
            n = len(entry_nested_list) * n
            for entry in entry_nested_list:
                if isinstance(entry, list):
                    assignWeightsToEntries(entry, weighted_entries, n)
                else:
                    weighted_entries.append((entry, 1 / float(n)))
            return weighted_entries

        kinetics, entry = self.estimateKineticsUsingRateRules(template)
        if entry:
            return [(entry, 1)], []  # Must be a rate rule
        else:
            # The template was estimated using an average or another node
            rules = []
            training = []

            lines = kinetics.comment.split('\n')

            # remove the Euclidean distance and family lines to help parser
            lines = [line for line in lines if not line.startswith('Euclid') and not line.startswith('family:')]

            # Discard the last line, unless it's the only line!
            # The last line is 'Estimated using ... for rate rule (originalTemplate)'
            # if from training reaction is in the first line append it to the end of the second line and skip the first line
            if 'Average of' not in kinetics.comment:
                if 'From training reaction' in lines[0]:
                    comment = lines[1]
                else:
                    comment = lines[0]
                if comment.startswith('Estimated using template'):
                    token_template_label = comment.split()[3][1:-1]
                    rule_entry, training_entry = self.retrieveOriginalEntry(token_template_label)
                    if training_entry:
                        training.append((rule_entry, training_entry, 1))  # Weight is 1
                    else:
                        rules.append((rule_entry, 1))
                else:
                    raise ValueError('Could not parse unexpected line found in kinetics comment: {}'.format(comment))
            else:
                comment = ' '.join(lines[:-1])
                # Clean up line for exec
                eval_comment_string = \
                    re.sub(r" \+ ", ",",                        # any remaining + signs
                    re.sub(r"Average of ", "",                  # average of averages
                    re.sub(r"Average of \[(?!Average)", "['",   # average of groups
                    re.sub(r"(\w|\))]", r"\1']",                # initial closing bracket
                    re.sub(r"(?<=[\w)]) \+ (?=Average)", "',",  # + sign between non-average and average
                    re.sub(r"(?<=]) \+ (?!Average)", ",'",      # + sign between average and non-average
                    re.sub(r"(?<!]) \+ (?!Average)", "','",     # + sign between non-averages
                    comment)))))))

                entry_nested_list = eval(eval_comment_string)

                weighted_entries = assignWeightsToEntries(entry_nested_list, [])

                rules = {}
                training = {}

                for token_template_label, weight in weighted_entries:
                    if 'From training reaction' in token_template_label:
                        token_template_label = token_template_label.split()[-1]
                    rule_entry, training_entry = self.retrieveOriginalEntry(token_template_label)
                    if training_entry:
                        if (rule_entry, training_entry) in training:
                            training[(rule_entry, training_entry)] += weight
                        else:
                            training[(rule_entry, training_entry)] = weight
                    else:
                        if rule_entry in rules:
                            rules[rule_entry] += weight
                        else:
                            rules[rule_entry] = weight
                # Each entry should now only appear once    
                training = [(k[0], k[1], v) for k, v in training.items()]
                rules = list(rules.items())

            return rules, training

    def extractSourceFromComments(self, reaction):
        """
        Returns the rate rule associated with the kinetics of a reaction by parsing the comments.
        Will return the template associated with the matched rate rule.
        Returns a tuple containing (Boolean_Is_Kinetics_From_Training_reaction, Source_Data)
        
        For a training reaction, the Source_Data returns::

            [Family_Label, Training_Reaction_Entry, Kinetics_In_Reverse?]
        
        For a reaction from rate rules, the Source_Data is a tuple containing::

            [Family_Label, {'template': originalTemplate,
                            'degeneracy': degeneracy,
                            'exact': boolean_exact?,
                            'rules': a list of (original rate rule entry, weight in average)
                            'training': a list of (original rate rule entry associated with training entry, original training entry, weight in average)}]


        where Exact is a boolean of whether the rate is an exact match, Template is
        the reaction template used, RateRules is a list of the rate rule entries containing
        the kinetics used, and TrainingReactions are ones that have created rules used in the estimate.
        """
        lines = reaction.kinetics.comment.split('\n')

        exact = False
        template = None
        rules = None
        training_entries = None
        degeneracy = 1

        regex = r"\[(.*)\]"  # only hit outermost brackets
        for line in lines:
            if line.startswith('Matched'):
                # Source of the kinetics is from training reaction
                training_reaction_index = int(line.split()[2])
                depository = self.getTrainingDepository()
                training_entry = depository.entries[training_reaction_index]
                # Perform sanity check that the training reaction's label matches that of the comments
                if training_entry.label not in line:
                    raise AssertionError('Reaction {0} uses kinetics from training reaction {1} '
                                         'but does not match the training reaction {1} from the '
                                         '{2} family.'.format(reaction, training_reaction_index, self.label))

                # Sometimes the matched kinetics could be in the reverse direction..... 
                if reaction.isIsomorphic(training_entry.item, eitherDirection=False):
                    reverse = False
                else:
                    reverse = True
                return True, [self.label, training_entry, reverse]

            elif line.startswith('Exact match'):
                exact = True
            elif line.startswith('Estimated'):
                pass
            elif line.startswith('Multiplied by'):
                degeneracy = float(line.split()[-1])

        # Extract the rate rule information 
        full_comment_string = reaction.kinetics.comment.replace('\n', ' ')

        # The rate rule string is right after the phrase 'for rate rule'
        rate_rule_string = full_comment_string.split("for rate rule", 1)[1].split()[0]

        if rate_rule_string[0] == '[':
            template_label = re.split(regex, rate_rule_string)[1]
        else:
            template_label = rate_rule_string  # if has the line 'From training reaction # for rate rule node1;node2'

        template = self.retrieveTemplate(template_label.split(';'))
        rules, training_entries = self.getSourcesForTemplate(template)

        if not template:
            raise ValueError('Could not extract kinetics source from comments for reaction {}.'.format(reaction))

        source_dict = {'template': template, 'degeneracy': degeneracy, 'exact': exact,
                       'rules': rules, 'training': training_entries}

        # Source of the kinetics is from rate rules
        return False, [self.label, source_dict]

    def getBackboneRoots(self):
        """
        Returns: the top level backbone node in a unimolecular family.
        """

        backbone_roots = [entry for entry in self.groups.top if entry in self.forwardTemplate.reactants]
        return backbone_roots

    def getEndRoots(self):
        """
        Returns: A list of top level end nodes in a unimolecular family
        """

        end_roots = [entry for entry in self.groups.top if entry not in self.forwardTemplate.reactants]
        return end_roots

    def getTopLevelGroups(self, root):
        """
        Returns a list of group nodes that are the highest in the tree starting at node "root".
        If "root" is a group node, then it will return a single-element list with "root".
        Otherwise, for every child of root, we descend until we find no nodes with logic
        nodes. We then return a list of all group nodes found along the way.
        """

        group_list = [root]
        all_groups = False

        while not all_groups:
            new_group_list = []
            for entry in group_list:
                if isinstance(entry.item, Group):
                    new_group_list.append(entry)
                else:
                    new_group_list.extend(entry.children)
            group_list = new_group_list
            all_groups = all([isinstance(entry.item, Group) for entry in group_list])

        return group_list


def informationGain(ks1, ks2):
    """
    calculates the information gain as the sum of the products of the standard deviations at each
    node and the number of reactions at that node
    """
    return len(ks1) * np.std(ks1) + len(ks2) * np.std(ks2)


def getObjectiveFunction(kinetics1, kinetics2, obj=informationGain, T=1000.0):
    """
    Returns the value of four potential objective functions to minimize
    Uncertainty = N1*std(Ln(k))_1 + N1*std(Ln(k))_1
    Mean difference: -abs(mean(Ln(k))_1-mean(Ln(k))_2)
    Error using mean: Err_1 + Err_2
    Split: abs(N1-N2)
    """
    ks1 = np.array([np.log(k.getRateCoefficient(T)) for k in kinetics1])
    ks2 = np.array([np.log(k.getRateCoefficient(T)) for k in kinetics2])
    N1 = len(ks1)

    return obj(ks1, ks2), N1 == 0


def _makeRule(rr):
    """
    function for parallelization of rule and uncertainty calculation
    Errors in Ln(k) at each reaction are treated as samples from a weighted normal distribution
    weights are inverse variance weights based on estimates of the error in Ln(k) for each individual reaction
    """
    recipe, rxns, Tref, fmax, label, ranks = rr
    n = len(rxns)
    for i, rxn in enumerate(rxns):
        rxn.rank = ranks[i]
    rxns = np.array(rxns)
    if n > 0:
        kin = ArrheniusBM().fitToReactions(rxns, recipe=recipe)
        if n == 1:
            kin.uncertainty = RateUncertainty(mu=0.0, var=(np.log(fmax) / 2.0) ** 2, N=1, Tref=Tref, correlation=label)
        else:
            dlnks = np.array([
                np.log(
                    ArrheniusBM().fitToReactions(rxns[list(set(range(len(rxns))) - {i})], recipe=recipe)
                    .toArrhenius(rxn.getEnthalpyOfReaction(Tref))
                    .getRateCoefficient(T=Tref) / rxn.getRateCoefficient(T=Tref)
                ) for i, rxn in enumerate(rxns)
            ])  # 1) fit to set of reactions without the current reaction (k)  2) compute log(kfit/kactual) at Tref
            varis = (np.array([rank_accuracy_map[rxn.rank].value_si for rxn in rxns]) / (2.0 * 8.314 * Tref)) ** 2
            # weighted average calculations
            ws = 1.0 / varis
            V1 = ws.sum()
            V2 = (ws ** 2).sum()
            mu = np.dot(ws, dlnks) / V1
            s = np.sqrt(np.dot(ws, (dlnks - mu) ** 2) / (V1 - V2 / V1))
            kin.uncertainty = RateUncertainty(mu=mu, var=s ** 2, N=n, Tref=Tref, correlation=label)
        return kin
    else:
        return None


def _spawnTreeProcess(family, templateRxnMap, obj, T, nprocs, depth, minSplitableEntryNum, minRxnsToSpawn):
    parent_conn, child_conn = mp.Pipe()
    name = list(templateRxnMap.keys())[0]
    p = mp.Process(target=_childMakeTreeNodes,
                   args=(family, child_conn, templateRxnMap, obj, T, nprocs,
                         depth, name, minSplitableEntryNum, minRxnsToSpawn))
    p.start()
    return parent_conn, p, name


def _childMakeTreeNodes(family, childConn, templateRxnMap, obj, T, nprocs, depth, name, minSplitableEntryNum,
                        minRxnsToSpawn):
    del_labels = []
    root_label = list(templateRxnMap.keys())[0]
    for label in family.groups.entries.keys():
        if label != root_label:
            del_labels.append(label)
    for label in del_labels:
        del family.groups.entries[label]

    family.groups.entries[root_label].parent = None

    family.makeTreeNodes(templateRxnMap=templateRxnMap, obj=obj, T=T, nprocs=nprocs, depth=depth + 1,
                         minSplitableEntryNum=minSplitableEntryNum, minRxnsToSpawn=minRxnsToSpawn)

    childConn.send(list(family.groups.entries.values()))
