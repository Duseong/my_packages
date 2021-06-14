from operator import add

def list_string_add(*args):
    addedstring = args[0]
    for arg in args[1:]:
        addedstring = list(map(add, addedstring, arg))

    return addedstring


