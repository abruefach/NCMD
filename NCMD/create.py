from ase.build.attach import attach_randomly_and_broadcast

def create_scene(atoms1, atoms2, distance, iterations):
    """
    Creates simulation scene of 2 different objects iteratively. Can be called
    multiple times to get a final scene with many structures.
    Parameters:
    Accepts:
    Returns:
    """
    for i in range(iterations):
        if i == 0:
            scene = attach_randomly_and_broadcast(atoms1, atoms2, distance)
        else:
            scene = attach_randomly_and_broadcast(scene, atoms2, distance)
    return scene