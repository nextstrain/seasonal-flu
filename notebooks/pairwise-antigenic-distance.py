import marimo

__generated_with = "0.13.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    from augur.utils import read_tree, read_node_data
    import json
    import Bio.SeqIO
    return Bio, json, read_node_data, read_tree


@app.cell
def _():
    tree_path = "builds/full-h1n1pdm/ha/tree.nwk"
    return (tree_path,)


@app.cell
def _(read_tree, tree_path):
    tree = read_tree(tree_path)
    return (tree,)


@app.cell
def _(tree):
    tree
    return


@app.cell
def _():
    subclades_path = "builds/full-h1n1pdm/ha/subclades.json"
    return (subclades_path,)


@app.cell
def _(read_node_data, subclades_path):
    subclades = read_node_data([subclades_path])["nodes"]
    return (subclades,)


@app.cell
def _(subclades, tree):
    # Find first node for each subclade.
    node_by_subclade = {}
    for node in tree.find_clades():
        subclade = subclades[node.name]["subclade"]
        if subclade not in node_by_subclade:
            node_by_subclade[subclade] = node.name
    return (node_by_subclade,)


@app.cell
def _(node_by_subclade):
    node_by_subclade
    return


@app.cell
def _(node_by_subclade):
    nodes = set(node_by_subclade.values())
    return (nodes,)


@app.cell
def _():
    ha1_sequences_path = "builds/full-h1n1pdm/ha/translations/HA1_withInternalNodes.fasta"
    return (ha1_sequences_path,)


@app.cell
def _(Bio, ha1_sequences_path, nodes):
    sequence_by_node = {
        record.name: str(record.seq)
        for record in Bio.SeqIO.parse(ha1_sequences_path, "fasta")
        if record.name in nodes
    }
    return (sequence_by_node,)


@app.cell
def _(sequence_by_node):
    sequence_by_node
    return


@app.cell
def _(node_by_subclade, sequence_by_node):
    sequence_by_subclade = {
        subclade: sequence_by_node[node]
        for subclade, node in node_by_subclade.items()
    }
    return (sequence_by_subclade,)


@app.cell
def _():
    titer_model_path = "builds/full-h1n1pdm/ha/titers-sub-model/cell_hi.json"
    return (titer_model_path,)


@app.cell
def _(json, titer_model_path):
    with open(titer_model_path, "r", encoding="utf-8") as fh:
        substitutions = json.load(fh)["substitution"]

    substitutions = {
        substitution.replace("HA1:", ""): weight
        for substitution, weight in substitutions.items()
    }
    return (substitutions,)


@app.cell
def _(substitutions):
    substitutions
    return


@app.function
def get_distance_between_clades(substitutions, clade_a_sequence, clade_b_sequence):
    distance = 0.0
    
    for i in range(len(clade_a_sequence)):
        if clade_a_sequence[i] != clade_b_sequence[i]:
            distance += substitutions.get(
                f"{clade_a_sequence[i]}{i + 1}{clade_b_sequence[i]}",
                0.0,
            )

    return distance


@app.cell
def _(sequence_by_subclade, substitutions):
    get_distance_between_clades(
        substitutions,
        sequence_by_subclade["C.1.9.3"],
        sequence_by_subclade["D.3.1"],
    )
    return


@app.cell
def _(sequence_by_subclade, substitutions):
    get_distance_between_clades(
        substitutions,
        sequence_by_subclade["D.3.1"],
        sequence_by_subclade["C.1.9.3"],
    )
    return


@app.cell
def _(sequence_by_subclade, substitutions):
    get_distance_between_clades(
        substitutions,
        sequence_by_subclade["C.1.9.3"],
        sequence_by_subclade["D"],
    )
    return


@app.cell
def _(sequence_by_subclade, substitutions):
    get_distance_between_clades(
        substitutions,
        sequence_by_subclade["C.1.9.3"],
        sequence_by_subclade["C.1.1"],
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
