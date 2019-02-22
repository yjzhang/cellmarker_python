# TODO: integrate cell ontology into the cellmarker db
# incorporate previous layers and stuff...

from owlready2 import get_ontology, default_world

default_world.set_backend(filename='cellmarker/data/cl.sqlite3')
onto = get_ontology('cl.owl').load()
# http://owlready.8326.n8.nabble.com/Accessing-class-by-its-name-in-owlready2-td457.html
namespace = onto.get_namespace("http://purl.obolibrary.org/obo/")
default_world.save()
classes = list(onto.classes())

for cl in classes:
    print(cl.label)
    print(cl.comment)
    print(cl.descendants())

# TODO: find a way to get parents of a node???
