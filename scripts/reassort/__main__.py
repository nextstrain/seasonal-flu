import argparse, sys, os, glob, json
import Bio
import Bio.SeqIO
import Bio.Phylo
from collections import defaultdict
from random import randint


def get_args():
    parser = argparse.ArgumentParser(
        description="Define reassortant clades",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--trees', nargs='+',
                        help="newick tree files")
    parser.add_argument('--distances', nargs='+',
                        help="distance metrics (for each tree)")
    parser.add_argument('--output',
                        help="json file containing clade definitions")

    args=parser.parse_args()
    return args


def write_node_data_json(tree, output_path, attr_map, missing_value = 0):
    strain_data = defaultdict(dict)
    for node in t.tree.find_clades():
        for attr, label in attr_map.items():
            strain_data[node.name][label] = getattr(node, attr, missing_value)
    with open(output_path, 'w') as fh:
        json.dump({"nodes": strain_data}, fh, indent=2, sort_keys = True)
    print("Saved node JSON with these attrs:", ",".join(attr_map.values()))
    print("Don't forget to add them to the auspice config JSON before running \"augur export\"")


class Trees():
    def __init__(self, treePaths, distancePaths):
        super().__init__()
        print("Trees class __init__")
        print("\tReading in these trees:", treePaths)
        self.tree = Bio.Phylo.read(treePaths[0], "newick") # <class 'Bio.Phylo.Newick.Tree'>
        self.set_parent_links_in_tree(self.tree)
        self.other_trees = [Bio.Phylo.read(p, "newick") for p in treePaths[1:]]
        self.strains_common_to_all_trees = self.get_common_strains_between_trees([self.tree, *self.other_trees])
        self.current_reassort_id = 0
        print("\tReading in distance files:", distancePaths)
        self.distances = []
        for p in distancePaths[1:]:
            with open(p) as fh:    
                self.distances.append(json.load(fh))

    @staticmethod
    def set_parent_links_in_tree(tree):
        for clade in tree.find_clades(order='level'):
            for child in clade:
                child.parent = clade
        tree.root.parent = False

    @staticmethod
    def get_common_strains_between_trees(trees):
        strains = [set([node.name for node in t.find_clades() if node.is_terminal()]) for t in trees]
        return set.intersection(*strains)

    def do_all_nodes_have_same_set_of_common_ancestors(self, node, comparison_trees):
        def equal_set_exists(a, listB):
            for b in listB:
                if a == b:
                    return True
            return False

        for n1 in node.find_clades(terminal=False, order='preorder'):
            set_of_strains = {n2.name for n2 in n1.get_terminals() if n2.name in self.strains_common_to_all_trees}
            # print("\tChecking set of strains from ", n1.name, "with ", len(set_of_strains), "strains")
            for terminal_sets_in_other_tree in self.terminal_sets:
                if not equal_set_exists(set_of_strains, terminal_sets_in_other_tree): 
                    return False # one failure means the question if false
        
        return True


    @staticmethod
    def mark_downstream(node, attr, value, include_node=True):
        if include_node:
            setattr(node, attr, value)
        for n in node.find_clades():
            setattr(n, attr, value)


    def store_terminal_sets(self):
        """
        for all trees, calculate & store sets of common tips across nodes.
        A convenient function for performance reasons.
        Does not consider terminal nodes.
        """
        self.terminal_sets = []
        for tree in self.other_trees:
            self.terminal_sets.append(
                [{n.name for n in node.get_terminals() if n.name in self.strains_common_to_all_trees} for node in tree.find_clades(terminal=False, order='preorder')]
            )
    
    def find_best_seeds(self):
        self.seed_nodes = []
        self.store_terminal_sets()
        seed_id = 1
        for node in self.tree.find_clades(terminal=False, order="preorder"):
            if hasattr(node, "seed_id"):
                continue
            if self.do_all_nodes_have_same_set_of_common_ancestors(node, self.other_trees):
                self.mark_downstream(node, "seed_id", seed_id)
                self.seed_nodes.append(node)
                print(node.name, node.count_terminals(), "has been assigned seed ID", seed_id)
                seed_id += 1


    def reset_pointer(self):
        """set the pointer to a new terminal (i.e. strain)
        that has not yet been assigned a reassort id.
        Reset the appropriate temp variables (variables prefixed with "current")
        """

        # things to do on every reset, no matter what
        self.current_reassort_id += 1
        self.current_false_moves = 0
        self.current_pruned_nodes = set([])
        self.current_pruned_taxa = set([])
        if not hasattr(self, "current_singleton_nodes"):
            self.current_singleton_nodes = set([])
        if len(self.current_singleton_nodes):
            for node in self.current_singleton_nodes:
                delattr(node, "reassort_id")
            self.current_singleton_nodes = set([])
        self.current_strain_set = set([])
        
        # try to find a new terminal node not part of a current reassort id clade
        for clade in self.tree.find_clades(order="postorder"):
            if not clade.is_terminal():
                continue
            if hasattr(clade, "reassort_id"):
                continue
            self.pointer = clade
            setattr(self.pointer, "reassort_id", self.current_reassort_id)
            print(f"\n\nPointer reset to a new seed: {self.pointer.name}. Reassort id: {self.current_reassort_id}")
            self.current_strain_set.add(self.pointer.name)
            return True

        print("\Can no longer reset pointer")
        return False

    def try_move_pointer(self):
        # if self.current_false_moves > 0:
        #     print("Not moving pointer -- current_false_moves limit exceeded")
        #     return False
        if not self.pointer.parent: # pointer is at the root
            print("Not moving pointer -- pointer already at root!")
            return False
        if hasattr(self.pointer.parent, "reassort_id"):
            print("Not moving pointer -- parent already assigned ID")
            return False
        self.pointer = self.pointer.parent
        return True


    def preorder_excluding_pruned(self, node):
        """preorder traversal excluding clades marked as pruned"""
        if node.is_terminal():
            yield node
        else:
            for child in node.clades:
                # if child.name in self.current_pruned_nodes:
                #     continue
                yield node
                yield from self.preorder_excluding_pruned(child)

    def _are_taxa_congruent(self, taxa_A):
        """MUST TAKE INTO ACCOUNT POLYTOMIES!!!"""

        taxa_A = taxa_A.intersection(self.strains_common_to_all_trees)
        taxa_A = taxa_A.difference(self.current_pruned_taxa)


        ## compare this to other trees! (this can be drastically sped up, BTW)
        # currently assume only 1 other tree to simplify things
        taxa_B = set([])
        for node in self.other_trees[0].common_ancestor(taxa_A).find_clades():
            if node.is_terminal:
                taxa_B.add(node.name)

        taxa_B = taxa_B.intersection(self.strains_common_to_all_trees)
        taxa_B = taxa_B.difference(self.current_pruned_taxa)

        # concordant if the taxa are the same!
        if len(taxa_A) == len(taxa_B):
            print("YES")
            return True
        print("NO ({} vs {})".format(len(taxa_A), len(taxa_B)))
        return False

    def _prune_child_clades(self, taxa_to_exclude):
        for clade in self.pointer.clades:
            if hasattr(clade, "reassort_id") and getattr(clade, "reassort_id") == self.current_reassort_id:
                pass
            else:
                self.current_pruned_nodes.add(clade)
        self.current_pruned_taxa = self.current_pruned_taxa.union(taxa_to_exclude)
        

    def try_add_children(self):
        """Given a pointer, one of the child clades will be the current reassort_id
        add the others 
        * IF they are terminal AND they preserve concordance between tree tip sets
        * IF they are a labelled clade AND they preserve concordance
        Probably involves some non-greedy algorithm for polytomies (not yet implemented -- currently greedy).
        Children added are also added to self.current_strain_set and reassort_id set accordingly
        """
        
        strains_to_consider = set([]) #self.current_strain_set.copy()
        nodes_to_maybe_add = []
        non_assigned_non_terminal_clade_is_child = False

        for clade in self.pointer.clades:
            if hasattr(clade, "reassort_id") and getattr(clade, "reassort_id") == self.current_reassort_id:
                pass
            elif clade.is_terminal():
                # doesn't matter if it's part of another reassort_id... we can steal it
                strains_to_consider.add(clade.name)
                nodes_to_maybe_add.append(clade)
            elif hasattr(clade, "reassort_id"):
                nodes_to_maybe_add.append(clade)
                for n in clade.find_clades():
                    nodes_to_maybe_add.append(n)
                    if n.is_terminal():
                        strains_to_consider.add(n.name)
            else:
                non_assigned_non_terminal_clade_is_child = True
        
        if non_assigned_non_terminal_clade_is_child:
            self.current_false_moves += 1
            print("not adding children as non_assigned_non_terminal_clade_is_child")
            self._prune_child_clades(taxa_to_exclude=strains_to_consider)
            self.current_singleton_nodes.add(self.pointer)
        elif len(nodes_to_maybe_add) and self._are_taxa_congruent(self.current_strain_set.union(strains_to_consider)):
            self.current_strain_set = self.current_strain_set.union(strains_to_consider)
            for node in nodes_to_maybe_add:
                setattr(node, "reassort_id", self.current_reassort_id)
            self.current_false_moves = 0
            self.current_singleton_nodes = set([])
            self.maybe_add_downstream_if_possible()

        else:
            self.current_false_moves += 1
            self._prune_child_clades(taxa_to_exclude=strains_to_consider)
            self.current_singleton_nodes.add(self.pointer)
        
        # in all cases, mark the pointer (sometimes it'll be in the current_singleton_nodes set)
        setattr(self.pointer, "reassort_id", self.current_reassort_id)


    def maybe_add_downstream_if_possible(self):

        # are any of the "pruned" taxa in the CA of tree 2?
        taxa_A = self.current_strain_set.intersection(self.strains_common_to_all_trees)
        terms_B = {node for node in self.other_trees[0].common_ancestor(taxa_A).get_terminals()}
        taxa_B = {node.name for node in terms_B}

        taxa_names_to_add = self.current_pruned_taxa.intersection(taxa_B)
        if not len(taxa_names_to_add):
            return
        
        print("ADDING BACK!", taxa_names_to_add)
        for node in self.pointer.get_terminals():
            if node.name in taxa_names_to_add:
                self.current_pruned_taxa.remove(node.name)
                setattr(node, "reassort_id", self.current_reassort_id)
                # linking nodes
                parent = node.parent
                while True:
                    if getattr(parent, "reassort_id", 0) == self.current_reassort_id:
                        break
                    setattr(parent, "reassort_id", self.current_reassort_id)
                    self.current_pruned_nodes.discard(parent) # no error if not present
                    parent = parent.parent

    def prune_singletons(self):
        for clade in self.tree.find_clades():
            if clade.is_terminal() and hasattr(clade, "reassort_id"):
                if getattr(clade, "reassort_id") != getattr(clade.parent, "reassort_id", "not_set"):
                    setattr(clade, "reassort_id", -1)



    def prune_twosomes(self):
        for clade in self.tree.find_clades():
            if len(clade.get_terminals()) <= 2 and hasattr(clade, "reassort_id"):
                if getattr(clade, "reassort_id") != getattr(clade.parent, "reassort_id", "not_set"):
                    setattr(clade, "reassort_id", -1)
                    for child in clade.clades:
                        setattr(child, "reassort_id", -1)

if __name__ == '__main__':
    args = get_args()

    t = Trees(args.trees, args.distances)
    debugging_break_counter = 0

    # Use phylogenetic concordance (strict) to identify the best seeds to start the algorithm
    t.find_best_seeds()

    write_node_data_json(t, args.output, {"seed_id": "seed_id"})

    # strain_data[node.name]["reassort_raw"] = getattr(node, "reassort_id", 0)
    # strain_data[node.name]["reassort"] = getattr(node, "reassort_id_no_gaps", 0)

    sys.exit(0)


    # Algorithm picks a leaf as a seed (and sets the pointer to it)
    # and we then try to greedily add related strains to the set,
    # according to if they have reassorted.
    while t.reset_pointer():

        debugging_break_counter += 1
        if debugging_break_counter>10:
            break

        while t.try_move_pointer():
            print("pointer moved up to {} ({} terminals)".format(t.pointer.name, len(t.pointer.get_terminals())))
            
            t.try_add_children()

            if t.current_false_moves > 0:
                print(f"Current false moves: {t.current_false_moves}")
                print(f"Excluded taxa: {t.current_pruned_taxa}")

            if t.current_false_moves > 5:
                print("TOO MANY FALSE MOVES")
                break


            # if t.is_pointer_concordant():
            #     t.mark_pointer()
            #     continue

            # if len(t.pointer.clades)>2:
            #     print("polytomy detected. Abandoning search")
            #     break

            # the children nodes of the pointer have broken concordance.
            # prune them off!

            # t.current_false_moves += 1
            # t.prune_children_from_pointer()
            # t.mark_clade_not_good()


    # set all nodes where 2 children & both are terminal & parent has different reassort_id
    # t.prune_twosomes()
    t.prune_singletons()

    
    ## trickery
    ids = {getattr(node, "reassort_id", 0) for node in t.tree.find_clades()}
    ids.discard(0)
    ids.discard(-1)
    print(ids)
    mapping = {old: idx+1 for idx, old in enumerate(list(ids))}
    mapping[-1] = -1
    print(mapping)
    for node in t.tree.find_clades():
        if hasattr(node, "reassort_id"):
            setattr(node, "reassort_id_no_gaps", mapping[getattr(node, "reassort_id")])

    # write output
    strain_data = defaultdict(dict)
    for node in t.tree.find_clades():
        strain_data[node.name]["reassort_raw"] = getattr(node, "reassort_id", 0)
        strain_data[node.name]["reassort"] = getattr(node, "reassort_id_no_gaps", 0)
    with open(args.output, 'w') as fh:
        json.dump({"nodes": strain_data}, fh, indent=2, sort_keys = True)

