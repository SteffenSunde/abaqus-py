def scalarFieldToVTK(odb, part_name, fields, step=-1, frame=-1, output_file = ''):
    """
    Directly generates a VTK Legacy file from the abaqus Output Database File.
    VTK file is appended point and element data, element centroid scalar fields.

    Input:
        odb         : Abaqus Odb object
        part_name   : Name of part to generate geometry from
        fields      : List of strings representing the fields to extract
        step        : (Optionally) Step at which fields are found
        frame       : (Optionally) Frame at which fields are found
        output_file : (Optional) The path name of the output VTK file

    Output:
        Generates an ASCII Legacy VTK file with Cauchy stress at element 
        centroids and nodal displacements.

    TODO:
        - Testing
        - Add vector fields
        - Add tensor fields

    """
    from abaqusConstants import CENTROID, SCALAR

    step = odb.steps.items()[step][1]
    part = odb.rootAssembly.instances[part_name + '-1']

    num_nodes = len(part.nodes)
    num_elements = len(part.elements)

    if output_file == '':
        vtk_path = part_name + '.vtk'
    else:
        vtk_path = output_file

    with open(vtk_path, 'w') as vtk_file:

        # For debugging
        print('Generating VTK file '+ vtk_path)

        # Write VTK header data
        vtk_file.write('# vtk DataFile Version 4.0')
        vtk_file.write('\nAbaqus part ' + part_name)
        vtk_file.write('\nASCII' )
        vtk_file.write('\nDATASET UNSTRUCTURED_GRID')

        # Append node coordinates to file
        print('Appending '+ str(num_nodes)+ ' node points to VTK file')
        vtk_file.write( '\n\nPOINTS ' + str(num_nodes) + ' double' )
        for i in range(num_nodes):
            node_position = part.nodes[i].coordinates
            vtk_file.write('\n')
            for j in range(3):
                vtk_file.write(str(node_position[j]) + ' ')

        #  Append cell connectivity to file
        cell_types = list()

        # Dictionary to translate Abaqus element types to corresponding VTK Legacy
        cell_type_dict = {'C3D8R': '12', 'C3D8': '12', 'C3D4': '10', 'C3D10': '24', 'CPS3': '5', 'CPE3': '5', 'CPS4R': '5', 'CPE4R': '5'}

        # Total number of data points in list of element connectivity
        num_cell_data_points = sum(len(element.connectivity) for element in part.elements)

        # Append elements to vtk file
        print('Appending {} elements points to VTK file'.format(str(num_elements)))
        vtk_file.write('\n\nCELLS ' + str(num_elements) + ' ' + str(num_cell_data_points+num_elements))
        for i in range(num_elements):
            cell_types.append(cell_type_dict[part.elements[i].type])
            connectivity = list(part.elements[i].connectivity)
            numNodesInElement = len(connectivity)
            vtk_file.write('\n')
            vtk_file.write(str(numNodesInElement)+' ')
            for pos in connectivity:
                vtk_file.write(str(pos-1) + ' ')

        # Append Cell types to file.
        vtk_file.write('\n\nCELL_TYPES '+ str(num_elements))
        for cellType in cell_types:
            vtk_file.write('\n'+cellType)

        # Append CENTROID cell data
        for field_name in fields:
            vtk_file.write('\n\nCELL_DATA '+ str(num_elements))
            vtk_file.write("\nSCALARS {} double".format(field_name.replace(" ", "_")))
            vtk_file.write('\nLOOKUP_TABLE default\n')

            field = step.frames[frame].fieldOutputs[field_name].getSubset(region=part, position=CENTROID, )
            print("Number of scalar field values in field {}: {}".format(field_name, str(len(field.values))))
            for scalar in field.values:            
                vtk_file.write(str(scalar.data) + '\n')
    pass