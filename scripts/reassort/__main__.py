import argparse, sys, os, glob, json
import Bio
import Bio.SeqIO
import Bio.Phylo
from collections import defaultdict
from random import randint
import numpy as np


def get_args():
    parser = argparse.ArgumentParser(
        description="Define reassortant clades",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--trees', nargs='+',
                        help="newick tree files")
    parser.add_argument('--mutations', nargs='+',
                        help="inferred mutations (for each tree)")
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
    print("Saved node JSON with these attrs:", ", ".join(attr_map.values()))
    print("Don't forget to add them to the auspice config JSON before running \"augur export\"")

def print_title(title):
    print("***************************************************************************")
    print(title)
    print("***************************************************************************")

def print_end_of_block():
    print("---------------------------------------------------------------------------\n\n")

class Trees():
    def __init__(self, params, treePaths, mutationPaths):
        super().__init__()
        print_title("Trees class __init__")
        self.params = type('params', (), {})()
        [setattr(self.params, k, v) for k, v in params.items()]
        print("\tReading in these trees:", treePaths)
        self.trees = [Bio.Phylo.read(p, "newick") for p in treePaths]
        self.name_to_node_maps = []
        for tree in self.trees:
            self.set_parent_links_in_tree(tree)
            self.name_to_node_maps.append(self.get_name_to_node_map(tree))
        self.strains_common_to_all_trees = self.get_common_strains_between_trees(self.trees)
        self.tree = self.trees[0] # for convenience
        self.current_reassort_id = 0
        self.reassort_id_strain_map = {}
        print("\tReading in mutation JSONs files:", mutationPaths)
        self.mutations = []
        for p in mutationPaths:
            with open(p) as fh:
                raw_data = json.load(fh)
                self.mutations.append(raw_data["nodes"])
        print_end_of_block()

    @staticmethod
    def set_parent_links_in_tree(tree):
        for clade in tree.find_clades(order='level'):
            for child in clade:
                child.parent = clade
        tree.root.parent = False
    
    @staticmethod
    def get_name_to_node_map(tree):
        data = {}
        for clade in tree.find_clades():
            data[clade.name] = clade
        return data

    @staticmethod
    def get_common_strains_between_trees(trees):
        strains = [set([node.name for node in t.find_clades() if node.is_terminal()]) for t in trees]
        return set.intersection(*strains)

    def do_all_nodes_have_same_set_of_common_ancestors(self, node, comparison_trees, monopyletic_tip_sets_per_tree):
        def is_a_in_b(a, listB):
            for b in listB:
                if a == b:
                    return True
            return False

        for n1 in node.find_clades(terminal=False, order='preorder'):
            set_of_strains = {n2.name for n2 in n1.get_terminals() if n2.name in self.strains_common_to_all_trees}
            # print("\tChecking set of strains from ", n1.name, "with ", len(set_of_strains), "strains")

            for monopyletic_tip_sets in monopyletic_tip_sets_per_tree:
                # check if our set of strains is present in some monophyly in the other tree
                if not is_a_in_b(set_of_strains, monopyletic_tip_sets): 
                    return False # one failure means the question is false
        
        return True


    @staticmethod
    def mark_downstream(node, attr, value, include_node=True):
        if include_node:
            setattr(node, attr, value)
        for n in node.find_clades():
            setattr(n, attr, value)

    def find_best_seeds(self):
        print_title("Finding best seeds")
        self.seed_nodes = []

        # for speed reasons, precalculate the list of each set of terminal nodes (done for every node in the tree)
        monopyletic_tip_sets_per_tree = []
        for tree in self.trees:
            monopyletic_tip_sets_per_tree.append(
                [{n.name for n in node.get_terminals() if n.name in self.strains_common_to_all_trees} for node in tree.find_clades(terminal=False, order='preorder')]
            )

        seed_id = 1
        for node in self.tree.find_clades(terminal=False, order="preorder"):
            if hasattr(node, "seed_id"):
                continue
            if self.do_all_nodes_have_same_set_of_common_ancestors(node, self.trees[1:], monopyletic_tip_sets_per_tree[1:]):
                self.mark_downstream(node, "seed_id", seed_id)
                self.seed_nodes.append(node)
                print("\t", node.name, node.count_terminals(), "has been assigned seed ID", seed_id)
                seed_id += 1
        print_end_of_block()

    def get_stats_on_best_seeds(self):
        """
        Get stats from the seeds (which are phylogenetically thought to be non-reassorting)
        which will be able to inform our parameter choices 
        """
        print_title("Statistics on the \"best\" seeds")
        
        mutation_distributions = [[] for _ in self.trees]
        num_tips_per_seed = []

        for seed_root in self.seed_nodes:
            tip_names = {n.name for n in seed_root.get_terminals() if n.name in self.strains_common_to_all_trees}
            num_tips_per_seed.append(len(tip_names))
            for idx, tree in enumerate(self.trees):
                nodes_in_path = self.minimal_nodes_to_form_path(tree, strains=tip_names)
                x = self.num_of_mutations_observed_along_path(self.mutations[idx], nodes_in_path)
                mutation_distributions[idx].append(x)


        msg = ""
        msg += "mean number of mutations across seed clades: {:.2f}\n".format(np.mean(mutation_distributions))
        msg += "per tree means: "
        for idx in range(0, len(self.trees)):
            msg += "{:.2f}, ".format(np.mean(mutation_distributions[idx]))
        msg += "\nRANGES. min: {} 25%: {:.1f} median: {:.1f} 75%: {:.1f} max: {}\n".format(np.min(mutation_distributions), np.percentile(mutation_distributions, 25), np.percentile(mutation_distributions, 50), np.percentile(mutation_distributions, 75), np.max(mutation_distributions))
        msg += "mean num tips per clade: {:.1f}\n".format(np.mean(num_tips_per_seed))
        msg += "mean mutations per additional tip: {:.2f}".format(np.mean(mutation_distributions)/np.mean(num_tips_per_seed))
        print(msg)
        print_end_of_block()



    def reset_pointer(self):
        """
        Set the pointer to a new position to begin (attempting to) add sister strains which have not reassorted
        This could be a node (if it's been identified as a seed node) OR
        a terminal node. Either way, it will not yet have been assigned a reassort id.
        
        SIDE EFFECTS:
        self.pointer                        points to new node (!)
        self.current_reassort_id            incremented by one
        node.reassort_id (in self.trees)

        RETURNS:
        bool
        """

        if self.current_reassort_id:
            # the (soon-to-be-previous) reassort_id may have been set on "extra" (higher) nodes
            # due to how we check for more nodes to possible add to it.
            # ensure that it's not set higher than the CA of the tips which have this ID
            n = self.tree.common_ancestor(self.reassort_id_strain_map[self.current_reassort_id])
            while True:
                if not n.parent:
                    break
                if getattr(n.parent, "reassort_id", False)==self.current_reassort_id:
                    delattr(n.parent, "reassort_id")
                n = n.parent


        self.current_reassort_id += 1
        self.reassort_id_strain_map[self.current_reassort_id] = set()

        for clade in self.tree.find_clades(order="preorder"):
            if hasattr(clade, "seed_id") and not hasattr(clade, "reassort_id"):
                self.pointer = clade
                setattr(self.pointer, "reassort_id", self.current_reassort_id)
                for n in self.pointer.find_clades():
                    setattr(n, "reassort_id", self.current_reassort_id)
                for n in self.pointer.get_terminals():
                    if n.name in self.strains_common_to_all_trees:
                      self.reassort_id_strain_map[self.current_reassort_id].add(n.name)
                print("pointer set to id#", self.current_reassort_id, "(seed id = ", getattr(clade, "seed_id"), ") num taxa:", len(self.reassort_id_strain_map[self.current_reassort_id]))
                return True

        for node in self.tree.get_terminals():
            if node.name not in self.strains_common_to_all_trees:
                continue
            if hasattr(node, "reassort_id"):
                continue
            self.pointer = node
            setattr(self.pointer, "reassort_id", self.current_reassort_id)
            self.reassort_id_strain_map[self.current_reassort_id].add(node.name)
            print("pointer set to id #", self.current_reassort_id, self.pointer.name)
            return True

        return False

    def move_pointer_up(self):
        if not self.pointer.parent: # pointer is at the root
            raise Exception("Not moving pointer -- pointer already at root!")
        if hasattr(self.pointer.parent, "reassort_id"):
            raise Exception("Not moving pointer -- parent already assigned ID")
        self.pointer = self.pointer.parent


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
    
    def reorder_id(self, new_attr):
        """
        given a tree with "reassort_id" set, create a new label (new_attr) which
        where 1 is the taxa set with the most tips, 2 is the next largest etc
        """
        id_counts = defaultdict(lambda: 0)
        for x in [getattr(tip, "reassort_id", -1) for tip in t.tree.get_terminals()]:
            if x > 0:
                id_counts[x] += 1
        ordered = sorted(id_counts.items(), key=lambda x: x[1], reverse=True)
        id_map = {id_count[0]: idx+1 for idx, id_count in enumerate(ordered)}
        for node in t.tree.find_clades():
            reassort_id = getattr(node, "reassort_id", 0)
            try:
                setattr(node, new_attr, id_map[reassort_id])
            except KeyError:
                setattr(node, new_attr, 0)


    @staticmethod
    def minimal_nodes_to_form_path(tree, common_ancestors=[], strains=[]):
        """
        which nodes are needed to form a path beetween "end-nodes"

        INPUTS:
        common_ancestors    list of sets of nodes which the CA should be used as "end-nodes"
        strain_sets         list of sets of nodes which should each be used as "end-nodes"
        note that nodes can be node objects or string names

        RETURNS:
        list of node objects
        """
        
        def as_names(taxa):
            """
            if taxa (in taxa_A or taxa_B) are node objects, they may have come from a different
            tree than the one provided. We care about names, so use node.names instead
            """
            ret = []
            # watch out -- a single Bio.Phylo.Newick.Clade object can be iterated through!
            assert(type(taxa).__name__ in ["list", "set"])
            for x in taxa:
                try:
                    ret.append(x.name)
                except AttributeError:
                    ret.append(x)
            return ret

        # note on get_path method: the target node is part of the returned list,
        # but the root of the tree is not.
        # because mutations (in augur) are defined as "those leading to the node"
        # this is not an issue

        paths = []
        for ca_nodes in common_ancestors:
            paths.append(tree.get_path(target=tree.common_ancestor(as_names(ca_nodes))))
        for end_nodes in as_names(strains):
            paths.append(tree.get_path(target = end_nodes))


        # names common to all sets (these aren't part of the true path, they're part of the
        # path between the root of the tree and the root of the path
        try:
            exclude = set.intersection(*[   {n.name for n in path} for path in paths])
        except:
            import pdb; pdb.set_trace()

        path_nodes = [] # the data this function returns
        node_names_added = set() # list of node names in path_nodes so we don't double up
        for path in paths:
            for node in path:
                if node.name not in exclude and node.name not in node_names_added:
                    path_nodes.append(node)
                    node_names_added.add(node.name)
        
        return path_nodes

    @staticmethod
    def num_of_mutations_observed_along_path(mutations, path):
        n_muts = 0
        for node in path:
            n_muts += len(mutations[node.name]["muts"])
        # print("\tn muts across path of len", len(path), "is", n_muts)
        return n_muts
        

    def does_clade_appear_non_reassortant(self, clade):
        """
        Given a clade to check, how many mutations are needed to add
        all members of the clade to the current reassort_id
        If this is over the parameter value "mutation_threshold_for_inclusion"
        then it's classified as a reassortant.

        Here we could add further phylogenetic checks in the future.
        """

        already_has_id = getattr(clade, "reassort_id", False)
        taxa_in_this_id = self.reassort_id_strain_map[self.current_reassort_id]

 

        if already_has_id:
            # we are attempting to incorporate in another clade with a reassort_id
            # note that the parent should not have a reassort_id
            if clade.parent and hasattr(clade.parent, "reassort_id"):
                raise Exception("does_clade_appear_non_reassortant parent has id")

            # we take all the strains in the already existing ID
            # there should be a better way 
            # we consider the path between the current taxa and all the ones in the already existing id
            # because the clade has already been assessed as non-reassorting
            # we instead consider the path between the current taxa set and the new clade
            strains = {n for n in clade.get_terminals() if getattr(n, "reassort_id", False) == already_has_id}
        else:
            # we take all termial nodes
            strains = {n for n in clade.get_terminals() if n.name in self.strains_common_to_all_trees}

        if already_has_id == 8:
            import pdb; pdb.set_trace()

        for idx, t in enumerate(self.trees):
            if already_has_id:
                nodes_in_path = self.minimal_nodes_to_form_path(t, common_ancestors=[taxa_in_this_id, strains])
            else:
                nodes_in_path = self.minimal_nodes_to_form_path(t, common_ancestors=[taxa_in_this_id], strains=strains)

            x = self.num_of_mutations_observed_along_path(self.mutations[idx], nodes_in_path)
            if x > self.params.mutation_threshold_for_inclusion:
                # print("\t\tEXCLUDED (", x, "mutations!)")
                return False

        # debugging / info messages
        if already_has_id:
            print("\tmerging id {} into id {} ({} strains)".format(already_has_id, self.current_reassort_id, len(strains)))
        else:
            msg = "\tadding {} strains to id {}".format(len(strains), self.current_reassort_id)
            if len(strains) == 1:
                msg += " ({})".format(next(iter(strains)).name)
            print(msg)

        return True


    def assign_id_to_node(self, node):
        setattr(node, "reassort_id", self.current_reassort_id)
        if node.is_terminal() and node.name in self.strains_common_to_all_trees:
            self.reassort_id_strain_map[self.current_reassort_id].add(node.name)

    def add_clade(self, clade):
        """
        For a given clade, set the reassort_id attr on it & it's descendants

        If the given clade has already been asigned a reassort_id, we replace it with the new one
        (note that this may skip some descendents, this is by design)
        If no reassort_id is set on the clade, we add it to all descendants.

        SIDE EFFECTS:
        downstream nodes                the "reassort_id" attr is set
        self.reassort_id_strain_map     downstream terminal node names added to set
                                        key removed if old id is overwritten
        """
        previous_id = getattr(clade, "reassort_id", False)

        # loop through all descendent children
        for node in clade.find_clades(order="preorder"):
            # if adding in a clade already assigned a reassort_id, skip descendents who may
            # not have that id
            if previous_id and getattr(node, "reassort_id", False) != previous_id:
                continue
            self.assign_id_to_node(node)

        # add the clade itself
        self.assign_id_to_node(clade)

        # if we've overwritten a previous_id, it no longer exists in the tree        
        if previous_id:
            del self.reassort_id_strain_map[previous_id]
            for node in self.tree.find_clades():
                if getattr(node, "reassort_id", False) == previous_id:
                    import pdb; pdb.set_trace()
                    raise Exception("should have overwritten id {} but it's still in the tree".format(previous_id))


    def recursively_add_sisters_to_seeds(self, current_skip_count=0):
        print("recursively_add_sisters_to_seeds (current_skip_count: {})".format(current_skip_count))
        ## Move pointer up one heirachy
        try:
            self.move_pointer_up()
        except Exception as err:
            print(err)
            return False
        
        ## For each "sister" clade, can we add it in according to our metrics?
        for clade in self.pointer.clades:
            if hasattr(clade, "reassort_id") and getattr(clade, "reassort_id") == self.current_reassort_id:
                # clade can be the previous pointer, which this catches
                continue  
            elif self.does_clade_appear_non_reassortant(clade):
                self.add_clade(clade)
            else:
                current_skip_count += 1 # want to see entire polytomy, not bail early

        # ensure that all strains in this reassort_id are linked pylogenetically
        self.assign_id_to_node(self.pointer)

        if current_skip_count > self.params.max_skip_count:
            print ("\tSkip count exceeded. Adding no more sisters")
            return False
    
        self.recursively_add_sisters_to_seeds(current_skip_count)


if __name__ == '__main__':
    args = get_args()
    params = {
        "max_skip_count": 4,
        "mutation_threshold_for_inclusion": 10
    }

    t = Trees(params, args.trees, args.mutations)

    # Use phylogenetic concordance (strict) to identify the best seeds to start the algorithm
    t.find_best_seeds()
    t.get_stats_on_best_seeds()

    # # Add strains to seeds!
    while t.reset_pointer():
        t.recursively_add_sisters_to_seeds()
        print_end_of_block()

    # clean up results for display
    t.prune_singletons()
    # t.prune_twosomes()
    t.reorder_id("reassort_reordered")

    write_node_data_json(t, args.output, {
        "seed_id": "seed_id",
        "reassort_id": "reassort_id",
        "reassort_reordered": "reassort_reordered"
    })

