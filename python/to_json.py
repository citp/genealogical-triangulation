import json
from sex import Sex

def to_json(population, labeled_nodes, related_nodes):
    all_nodes = list(population.members)
    all_nodes.sort(key = lambda node: node._id)
    generation_map = population.node_to_generation
    json_nodes = list()
    for node in all_nodes:
        json_node = dict()
        json_node["id"] = node._id
        json_node["generation"] = generation_map[node]
        if node.sex == Sex.Male:
            json_node["sex"] = "Male"
        else:
            json_node["sex"] = "Female"
        if node.father is not None:
            json_node["father"] = node.father._id
        if node.mother is not None:
            json_node["mother"] = node.mother._id
        if node.suspected_father is not None:
            json_node["suspected_father"] = node.suspected_father._id
        if node.suspected_mother is not None:
            json_node["suspected_mother"] = node.suspected_mother._id
        if node.twin is not None:
            json_node["twin"] = node.twin._id
        json_nodes.append(json_node)
    json_labeled_nodes = [node._id for node in labeled_nodes]
    labeled_set = set(json_labeled_nodes)
    
    json_related = list()
    for a, b in related_nodes:
        if a in labeled_set:
            to_append = {"labeled_node": a._id,
                         "unlabeled_node": b._id}
        else:
            to_append = {"labeled_node": b._id,
                         "unlabeled_node": a._id}
        json_related.append(to_append)

    population_json = {"nodes": json_nodes,
                       "related": json_related,
                       "labeled": json_labeled_nodes}
    return json.dumps(population_json)
    
