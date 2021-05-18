from ase.build.attach import attach_randomly_and_broadcast

def create_scene(atoms1, atoms2, distance, iterations):
    """
    Creates simulation scene of 2 different objects iteratively. Can be called
    multiple times to get a final scene with many structures.
    Parameters:
    Accepts:
    atoms1        (Atoms Object) First atoms object. Can be an individual item or scene.
    atoms2        (Atoms Object) Second atoms object. Can be an individual item or scene.
                         Will iteratively attach this object to the first based on distance
                         and iterations inputs
    distance      (int) Minimum distance that atoms1 and atoms2 objects will be placed.
    iterations    (int) Number of iterations to perform
    Returns:
    scene         (Atoms Object) Final scene
    """
    for i in range(iterations):
        if i == 0:
            scene = attach_randomly_and_broadcast(atoms1, atoms2, distance)
        else:
            scene = attach_randomly_and_broadcast(scene, atoms2, distance)
    return scene