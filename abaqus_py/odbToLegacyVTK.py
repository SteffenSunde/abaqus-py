def scalar_field_to_vtk(odb, part_name, fields, step=-1, frame=-1, output_file = ''):
    """
    Directly generates a VTK Legacy file from the abaqus Output Database File.
    VTK file is appended point and element data, element centroid scalar fields.

    Input:
        odb         : Abaqus Odb object
        part_name   : Name of part to generate geometry from
        fields      : List of strings representing the fields to extract
        step        : (Optional) Step at which fields are found
        frame       : (Optional) Frame at which fields are found
        output_file : (Optional) The path name of the output VTK file

    Output:
        Generates an ASCII Legacy VTK file with Cauchy stress at element 
        centroids and nodal displacements.

    TODO:
        - Testing
        - Add vector fields
        - Add tensor fields
        - Add field data for subsets (e.g. 0 for the remaining points)

    """
    from abaqusConstants import CENTROID, SCALAR

    step_handle = odb.steps.items()[step][1]
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
        cell_type_dict = {'C3D8R': '12', 'C3D8': '12', 'C3D4': '10', 'C3D10': '24', 'CPS3': '5', 'CPE3': '5', 'CPS4R': '5', 'CPE4R': '5', 'CPE8R': '23', 'CPS8R': '23'}

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

            field = step_handle.frames[frame].fieldOutputs[field_name].getSubset(region=part, position=CENTROID, )
            print("Number of scalar field values in field {}: {}".format(field_name, str(len(field.values))))
            for scalar in field.values:            
                vtk_file.write(str(scalar.data) + '\n')
    return


def mesh_to_legacy_vtk(odb, part_name, output_file = '', step='', frame=-1):
    """
    Directly generates a VTK Legacy file from the abaqus Output Database File.
    VTK file is appended point and element data, element centroid stress tensor
    and nodal displacement vectors.

        Input:
            odb         : Abaqus Odb object
            part_name   : Name of part to generate geometry from (instance name)
            output_file : (Optional) The path name of the output VTK file
            step        : (Optional) Step at which fields are found
            frame       : (Optional) Frame at which fields are found

        Output:
            Generates an ASCII Legacy VTK file with Cauchy stress at element 
            centroids and nodal displacements.

        TODO:
            -
    """
    try:
        # from odbAccess import openOdb
        from abaqusConstants import CENTROID, NODAL
        part = odb.rootAssembly.instances[part_name.upper()]
    except ImportError as ex:
        print("Error importing from Abaqus library: " + str(ex))
    except KeyError as ex:
        instances = ""
        for key, _ in odb.rootAssembly.instances.items():
            instances += " \'"+ key + "\',"
        raise IOError("Could not find instance in odb-file: " + str(ex) +". Instance(s) found: " + instances)

    num_nodes = len(part.nodes)
    num_elements = len(part.elements)

    if output_file == '':
        vtk_path = part_name + '.vtk'
    else:
        vtk_path = output_file

    with open(vtk_path, 'w') as vtk_file:

        if step == '':
            step = odb.steps.items()[0][0]

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
        cell_type_dict = {'C3D8R': '12', 'C3D8': '12', 'C3D4': '10', 'C3D10': '24', 'CPS3': '5', 'CPE3': '5', 'CPS4R': '5', 'CPE4R': '5', 'CPE8R': '23', 'CPS8R': '23'}

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

        findley_csv_file = 'C:/temp/conf/QuadQuad/max_findley.csv'
        frame = odb.steps[step].frames[frame]
        field = frame.fieldOutputs['S'].getSubset(region=part, position=CENTROID)
        vtk_file.write('\n\nCELL_DATA '+ str(len(field.values)))
        csv_to_legacy_vtk(vtk_file, findley_csv_file)
        cellDataToLegacyVTK(vtk_file, field, frame)
        cellDataToLegacyVTK(vtk_file, field, frame, invariant="MISES")
        
        field = frame.fieldOutputs['LE'].getSubset(region=part, position=CENTROID)
        cellDataToLegacyVTK(vtk_file, field, frame)

        # Write nodal stress tensors
        # Abaqus does not record nodal stresses by default, either activate it,
        # or interpolate average the stresses from integration points?
        vtk_file.write('\n\nPOINT_DATA '+ str(num_nodes))
        #vtk_file.write('\nTENSORS NodalCauchyStress double')
        #vtk_file.write('\nLOOKUP_TABLE default')
        field_values = frame.fieldOutputs['S'].getSubset(region=part,
                                                        position=NODAL).values

        # Append global node stress tensors. IF they exsists. Abaqus
        # does not record global node stresses by default, but can be
        # activated using These are averaged from the
        # surrounding elements.
        if len(field_values) == num_nodes:
                for tensor in field_values:
                    vtk_file.write('\n' + str(tensor.data[0]) + ' ') # Sxx
                    vtk_file.write(str(tensor.data[3]) + ' ') # Sxy
                    vtk_file.write(str(tensor.data[4])) # Sxz
                    vtk_file.write('\n' + str(tensor.data[3]) + ' ') # Syx
                    vtk_file.write(str(tensor.data[1]) + ' ') # Syy
                    vtk_file.write(str(tensor.data[5])) # Syz
                    vtk_file.write('\n' + str(tensor.data[4]) + ' ') # Szx
                    vtk_file.write(str(tensor.data[5]) + ' ') # Szy
                    vtk_file.write(str(tensor.data[2]) + '\n') # Szz

        # Append displacement vectors at nodes.
        vtk_file.write('\nVectors NodalDisplacements float')
        displacement_field = frame.fieldOutputs['U'].getSubset(region=part,
                                                        position=NODAL).values
        for vector in displacement_field:
            if len(vector.data) == 3:
                vtk_file.write('\n'+ str(vector.data[0]) + ' ' + str(vector.data[1]) + ' ' + str(vector.data[2]))
            elif len(vector.data) == 2:
                vtk_file.write('\n'+ str(vector.data[0]) + ' ' + str(vector.data[1]) + ' 0.0')

    return 


def cellDataToLegacyVTK(output_handle, field_handle, frame_handle, invariant=''):
    """
    Outputs the field to the given output handle (VTK Legacy). 
    Be aware of the order of fields applied!

    TODO:
    - Testing
    - Add TRESCA invariant
    - Add MAX_PRINCIPAL invariant
    - Add MIN_PRINCIPAL invariant
    """
    field_name = field_handle.name.replace(" ", "_")
    field_type = str(field_handle.type)[0:6]        # SCALAR | VECTOR | TENSOR

    if field_type == "TENSOR" and invariant == "":
        output_handle.write("\nTENSORS {} double".format(field_name))
        print("Appending cell tensor field {}:".format(str(len(field_handle.values))))
        if len(field_handle.componentLabels) == 6: # 3D tensor
            for tensor in field_handle.values: 
                output_handle.write('\n' + str(tensor.data[0]) + ' ') # Sxx
                output_handle.write(str(tensor.data[3]) + ' ') # Sxy
                output_handle.write(str(tensor.data[4])) # Sxz
                output_handle.write('\n' + str(tensor.data[3]) + ' ') # Syx
                output_handle.write(str(tensor.data[1]) + ' ') # Syy
                output_handle.write(str(tensor.data[5])) # Syz
                output_handle.write('\n' + str(tensor.data[4]) + ' ') # Szx
                output_handle.write(str(tensor.data[5]) + ' ') # Szy
                output_handle.write(str(tensor.data[2]) + '\n') # Szz
        elif len(field_handle.componentLabels) == 4: # 2D tensor
            for tensor in field_handle.values:
                output_handle.write('\n' + str(tensor.data[0]) + ' ') # Sxx
                output_handle.write(str(tensor.data[3]) + ' ') # Sxy
                output_handle.write('0.0') # Sxz
                output_handle.write('\n' + str(tensor.data[3]) + ' ') # Syx
                output_handle.write(str(tensor.data[1]) + ' ') # Syy
                output_handle.write('0.0') # Syz
                output_handle.write('\n' + '0.0 ') # Szx
                output_handle.write('0.0 ') # Szy
                output_handle.write(str(tensor.data[2]) + '\n') # Szz
    elif field_type == "TENSOR" and invariant == "MISES":
            output_handle.write("\nSCALARS MISES double")
            output_handle.write('\nLOOKUP_TABLE default\n')
            print('Appending Von Mises field values:')
            for value in field_handle.values:            
                output_handle.write(str(value.mises) + '\n')
    elif field_type == "VECTOR": # Are there really any interesting cell vectors?
        output_handle.write("\nVectors {} float".format(field_name))
        for vector in field_handle.values:
            if len(vector.data) == 3:
                output_handle.write('\n'+ str(vector.data[0]) + ' ' + str(vector.data[1]) + ' ' + str(vector.data[2]))
            elif len(vector.data) == 2:
                output_handle.write('\n'+ str(vector.data[0]) + ' ' + str(vector.data[1]) + ' 0.0')
    elif field_type == "SCALAR":
        output_handle.write("\nSCALARS {} double".format(field_name))
        output_handle.write('\nLOOKUP_TABLE default\n')

        print("Appending cell scalar field {}:".format(field_name))
        for value in field_handle.values:            
            output_handle.write(str(value.mises) + '\n')

    return


def nodalDataToLegacyVTK(output_handle, field_handle, frame_handle, invariant=''):
    field_name = field_handle.name.replace(" ", "_")
    field_type = str(field_handle.type)[0:6]        # SCALAR | VECTOR | TENSOR
    
    if field_type == "VECTOR":
        # Append displacement vectors at nodes.
        output_handle.write('\nVectors {} float'.format(field_name))

        for vector in field_handle.values:
            if len(vector.data) == 3:
                output_handle.write('\n'+ str(vector.data[0]) + ' ' + str(vector.data[1]) + ' ' + str(vector.data[2]))
            elif len(vector.data) == 2:
                output_handle.write('\n'+ str(vector.data[0]) + ' ' + str(vector.data[1]) + ' 0.0')
    return


def csv_to_legacy_vtk(output_handle, csv_file_name):


    """
    Translates a csv field to vtk legacy cell data..

    TODO:
    - Support for partial fields
    """

    import csv

    with open(csv_file_name, 'r') as csv_file:
        reader = csv.reader(csv_file)
        num_rows = sum([1 for row in reader])

        if num_rows < 0:
            raise Exception("File empty!")
        
        csv_file.seek(0)
        #header = next(reader)
        field_name = 'Findley'
        #print("Header: {}".format(header))
        
        print("Number of rows: {}".format(num_rows))
    	
        output_handle.write('\n\nCELL_DATA {}'.format(num_rows))
        output_handle.write('\nSCALARS {} double'.format(field_name))
        output_handle.write('\nLOOKUP_TABLE default\n')
    
        for row in reader:
            output_handle.write('{}\n'.format(row[1]))
    return


