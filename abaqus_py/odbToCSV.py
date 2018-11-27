def field_to_csv(odb_handle, elset_handle, field_name, step_number=-1, frame_number=-1):
    """
    Outputs the given field at given frame number in odb file.

    TODO:
    - Make similar for nodal quantities
    - Test
    """
    from abaqusConstants import CENTROID
    import csv

    frame_handle = odb_handle.steps.items()[step_number][1].frames[frame_number]
    field = frame_handle.fieldOutputs[field_name].getSubset(region=elset_handle, position=CENTROID)
    file_name = 'CENTROID-{}-{}-{}-{}.csv'.format(elset_handle.name, field.name, step_number, frame_number)
    print("Exporting {}".format(file_name))
    with open(file_name, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')
        csv_writer.writerow(('Step time:{}'.format(frame_handle.frameValue),))
        csv_writer.writerow(('Element_id',) + field.componentLabels)
        for point in field.values:
            csv_writer.writerow((point.elementLabel,) + tuple(point.data))


def field_history_to_csv(odb_handle, elset_handle, field_name, steps=[-1], frequency=1):
    """
    Creates one csv file per frame distributed with given frequency in the given step.

    TODO:
    - Test
    - Perhaps create a folder?
    """

    for step in steps:
        step_handle = odb_handle.steps.items()[step][1]
        frames = distributed_selection(range(len(step_handle.frames)), frequency)
        #frames = [step.frames[-1]]
        for frame in frames:
            field_to_csv(odb_handle, elset_handle, field_name, step, frame)


def distributed_selection(lst, num):
    """
    Selects evenly elements from the given list with the given number of elements.
    Ensures that the last point are included. The first point is not necessary due
    to being part of the previous time step.

    For two points, the mid and last point should be selected...
    
    TODO:
    - Test
    """

    if num >= len(lst):
        return lst

    selection = []

    if num > 1:
        selection = [lst[1]] # Yes, not the first element.
        if num < len(lst):
            #q, _ = divmod(len(lst), num-1)
            q = len(lst)//(num-1)
            selection += [lst[i] for i in range(q, len(lst)-1, q)]

    selection.append(lst[-1])     

    return selection


def csv_to_legacy_vtk(csv_file_name, vtk_file_name):
    """
    Translates a csv field to vtk legacy cell data..

    TODO:
    - Support for partial fields
    """

    import csv

    with open(csv_file_name, 'r') as csv_file:
        reader = csv.reader(csv_file)
        num_rows = sum([1 for row in reader]) - 1

        if num_rows < 0:
            raise Exception("File empty!")
        
        csv_file.seek(0)
        header = next(reader)
        field_name = header[2]
        print("Header: {}".format(header))
        
        print("Number of rows: {}".format(num_rows))
    	
        with open(vtk_file_name, 'w+') as vtk_file:
            vtk_file.write('\n\nCELL_DATA {}'.format(num_rows))
            vtk_file.write('\nSCALARS {}'.format(field_name))
            vtk_file.write('\nLOOKUP_TABLE default\n')
        
            for row in reader:
                vtk_file.write('{}\n'.format(row[1]))
    return


def element_stress_history_to_csv(odb, region, output_file='', steps='', frame=-1):
    """
    Extracts centroid stress history for an element set and stores in comma separated file.

    TODO:
    - Read field once for each time step? Might be more efficient.
    """

    try:
        from abaqusConstants import CENTROID
        element_set = odb.rootAssembly.elementSets[region.upper()]
    except ImportError as ex:
        raise ImportError('Failed to import abaqus modules!')
    except KeyError as ex:
        raise IOError('Could not find element set in odb-file: ' + str(ex) + '. Element set(s) found: ' + str([el[0] for el in odb.rootAssembly.elementSets.items()]))

    if output_file == '':
        csv_path = region + '.csv'
    else:
        csv_path = output_file

    if steps == '':
        step_keys = odb.steps.keys()
    else:
        step_keys = [odb.steps.keys()[i] for i in steps]

    # For debugging
    print('Found {} elements in region {}'.format(len(element_set.elements[0]), region))
    print('Number of history points: {}'.format(len(step_keys)))
    print('Generating CSV file '+ csv_path)

    if element_set.elements[0][0].type[1] == '3':
        dimension = '3D'
    else:
        dimension = '2D'

    with open(csv_path, 'w') as csv_file:
        csv_file.write('Centroid {} stress tensor history\n'.format(dimension))

        for element in element_set.elements[0]:
            row = str(element.label) + ','
            for step in step_keys:
                field = odb.steps[step].frames[frame].fieldOutputs['S'].getSubset(region=element, position=CENTROID) # Extrapolated and averaged?
                for component in field.values[0].data:
                    row = row + str(component) + ','
            row = row + '\n'
            csv_file.write(row)

    return


def element_data_to_csv(odb_handle, part_name, step_number=-1, frame_number=-1):
    """
    Extracts data for fatigue post-processing routines.
    Each row represents an element.

    Row data contains:
     - Element label
     - Element centroid position
     - Stress tensor
     - Strain tensor
     - Strain energy density

    TODO:
    - Fix step time in output file
    - Make similar for nodal quantities
    - Test 
    """
    from abaqusConstants import CENTROID, SCALAR
    import csv

    instance_handle = odb_handle.rootAssembly.instances[part_name.upper()]
    step_handle =  odb_handle.steps.items()[step_number][1]
    frame_handle = step_handle.frames[frame_number]
    file_name = 'CENTROID-{}-{}-{}.csv'.format(part_name, step_number, frame_number)
    stress_field = frame_handle.fieldOutputs['S'].getSubset(region=instance_handle, position=CENTROID)
    strain_field = frame_handle.fieldOutputs['LE'].getSubset(region=instance_handle, position=CENTROID)  # "LE" is outputted for NLGEOM=ON
    strain_energy_density_field = frame_handle.fieldOutputs['SENER'].getSubset(region=instance_handle, position=CENTROID)
    total_time = step_handle.totalTime + frame_handle.frameValue

    print("Exporting {}".format(file_name))
    with open(file_name, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')
        csv_writer.writerow(['Step time:{}'.format(total_time)])
        # Should label row start with comma?
        csv_writer.writerow(['Element_id','x', 'y', 'z'] + list(stress_field.componentLabels) + list(strain_field.componentLabels) + [strain_energy_density_field.name])
        for i in range(len(stress_field.values)):
            element_label = stress_field.values[i].elementLabel
            position = centroid_position(instance_handle, element_label)
            csv_writer.writerow([element_label] + position + list(stress_field.values[i].data) + list(strain_field.values[i].data) + [strain_energy_density_field.values[i].data])


def centroid_position(instance_handle, element_label):
    """
    Note! The element (and node) labels in Abaqus starts with 1. 
    Be careful not to mix indices and labels.
    """
    node_ids = instance_handle.elements[element_label - 1].connectivity

    coordinate_sum = [0, 0, 0]
    try:
        # Warning: Labels are 1-indexed
        # Label id is not necessarily same as node indices (See Abaqus guide 34.17.1)
        for node in [instance_handle.nodes[n - 1] for n in node_ids]:  # Sure this is safe?
            coordinate_sum += node.coordinates
    except Exception as err:
        print("Error reading element {} having nodes {}: {}".format(element_label, str(node_ids), err))

    return [coord / len(node_ids) for coord in coordinate_sum]


def element_at_pos(instance_handle, x, y):
    pass


def nodal_data(field_name, instance_handle, frame_handle):
    # Loop through elements
        # Loop through connectivity
    # For each node:
     # Store counter for number of elements
     # Store stress contribution for each element
     # See Abaqus guide "54.3.2 Selecting field output variables"
     # and "42.6.1 Understanding how results are computed"

     # TODO: 
     #  - Check values
     #  - Apply averaging scheme? See guide 42.6.2 Understanding result value averaging

    # Average
    import numpy as np
    from abaqusConstants import ELEMENT_NODAL
    
    # Depends on the field being ordered with ascending node label? No
    field = frame_handle.fieldOutputs[field_name].getSubset(region=instance_handle, position=ELEMENT_NODAL)

    num_nodes = len(instance_handle.nodes)
    num_components = max(len(field.componentLabels), 1)

    result = np.zeros((num_nodes, num_components))  # Initialized to zero?
    connected_elements = np.zeros((num_nodes, 1))  # Initialized to zero?

    for element_node in field.values:
        node_label = element_node.nodeLabel  
        result[node_label - 1, :] += element_node.data
        connected_elements[node_label - 1] += 1

    if field.componentLabels == ():
        labels = [field.name]
    else:
        labels = list(field.componentLabels)
    
    # TODO Is this correct? Important
    for i in range(num_nodes):
        result[i] /= connected_elements[i]

    return result, labels


def nodal_data_to_csv(odb_handle, part_name, step_number=-1, frame_number=-1):
    """
    Extracts fatigue data at nodes. Nodal values are obtained by extrapolating
    element integration point results and averaging each elements contribution.

    TODO: Check SENER output!

    """
    from abaqusConstants import SCALAR
    import csv

    instance_handle = odb_handle.rootAssembly.instances[part_name.upper()]
    step_handle =  odb_handle.steps.items()[step_number][1]
    frame_handle = step_handle.frames[frame_number]
    file_name = 'NODAL-{}-{}-{}.csv'.format(part_name, step_number, frame_number)
    stress_field, stress_labels = nodal_data('S', instance_handle, frame_handle)
    strain_field, strain_labels = nodal_data('LE', instance_handle, frame_handle)
    sener_field, sener_label = nodal_data('SENER', instance_handle, frame_handle)
    total_time = step_handle.totalTime + frame_handle.frameValue

    print("Exporting {}".format(file_name))
    with open(file_name, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')
        csv_writer.writerow(['Step time:{}'.format(total_time)])
        csv_writer.writerow(['Node_label', 'x', 'y', 'z'] + stress_labels + strain_labels + sener_label)
        for node in instance_handle.nodes:
            node_label = node.label
            position = instance_handle.nodes[node_label - 1].coordinates
            csv_writer.writerow([node.label] + list(position) + list(stress_field[node_label - 1,:]) + list(strain_field[node_label - 1, :]) + list(sener_field[node_label - 1]))


###################################################

# TODO 21.11.2018: Implement these functions to work for
# a subset of nodes or elements given by name (Geometric set preferably!)
# Probably go through the elementSet.

def nodal_data_section(field_name, odb_handle, set_name, frame_handle, csys_name=""):
    # Loop through elements
        # Loop through connectivity
    # For each node:
     # Store counter for number of elements
     # Store stress contribution for each element
     # See Abaqus guide "54.3.2 Selecting field output variables"
     # and "42.6.1 Understanding how results are computed"

     # TODO: 
     #  - Check values
     #  - Apply averaging scheme? See guide 42.6.2 Understanding result value averaging

    # Average
    import numpy as np
    from abaqusConstants import ELEMENT_NODAL
    
    # Depends on the field being ordered with ascending node label? No
    set_handle = odb_handle.rootAssembly.nodeSets[set_name.upper()]
    if csys_name == "":
        field = frame_handle.fieldOutputs[field_name].getSubset(region=set_handle, position=ELEMENT_NODAL)
    else:
        csys_handle = odb_handle.rootAssembly.datumCsyses[csys_name.upper()]
        field = frame_handle.fieldOutputs[field_name].getSubset(region=set_handle, position=ELEMENT_NODAL).getTransformedField(datumCsys=csys_handle)

    num_components = max(len(field.componentLabels), 1)
    result_dict = {}
    connected_dict = dict([(node.nodeLabel, 0) for node in field.values])

    # This doesnt work because for node set 
    for element_node in field.values: #  len(field.values) is not len(nodes) ??
        node_label = element_node.nodeLabel  
        # result[i, :] += element_node.data
        # connected_elements[i] += 1
        result_dict[node_label] = element_node.data
        connected_dict[node_label] += 1
        #result[node_label - 1, :] += element_node.data  # How to name these with labels of arbitrary integer id? TODO
        #connected_elements[node_label - 1] += 1

    if field.componentLabels == ():
        labels = [field.name]
    else:
        labels = list(field.componentLabels)
    
    # Remove the nodeLabel at last position?
    # Will the node ids always be correct over different fields? TODO
    result = np.ndarray((len(result_dict), num_components))
    for i, node in enumerate(set_handle.nodes[0]):
        node_label = node.label
        result[i, :] = np.asarray(result_dict[node_label]) / connected_dict[node_label]
        # result[i, -1] = node_label  # TODO

    return result, labels


def nodal_data_to_csv_section(odb_handle, set_name, step_number=-1, frame_number=-1, csys_name=""):
    """
    TODO: Check SENER output!
    TODO: Fix for arbitrary node set!

    """
    from abaqusConstants import SCALAR
    import csv
    import numpy as np

    if csys_name != "":
        csys_handle = odb_handle.rootAssembly.datumCsyses[csys_name.upper()]
        transformation =    np.outer(csys_handle.xAxis, [1, 0, 0]) + \
                            np.outer(csys_handle.yAxis, [0, 1, 0]) + \
                            np.outer(csys_handle.zAxis, [0, 0, 1])
        transformation = transformation.transpose()
        translation = csys_handle.origin

    # nodal_data_section(field_name, odb_handle, set_name, frame_handle):
    #instance_handle = odb_handle.rootAssembly.instances[part_name.upper()]
    set_handle = odb_handle.rootAssembly.nodeSets[set_name.upper()]
    step_handle =  odb_handle.steps.items()[step_number][1]
    frame_handle = step_handle.frames[frame_number]
    file_name = 'NODAL-{}-{}-{}.csv'.format(set_name, step_number, frame_number)
    stress_field, stress_labels = nodal_data_section('S', odb_handle, set_name, frame_handle, csys_name)
    strain_field, strain_labels = nodal_data_section('LE', odb_handle, set_name, frame_handle, csys_name)
    sener_field, sener_label = nodal_data_section('SENER', odb_handle, set_name, frame_handle)
    total_time = step_handle.totalTime + frame_handle.frameValue

    # For coordinates:
    # Find translation?
    # Rotate?

    print("Exporting {}".format(file_name))
    with open(file_name, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')
        csv_writer.writerow(['Step time:{}'.format(total_time)])
        csv_writer.writerow(['Node_label', 'x', 'y', 'z'] + stress_labels + strain_labels + sener_label)
        for i, node in enumerate(set_handle.nodes[0]):
            node_label = node.label
            #position = instance_handle.nodes[node_label - 1].coordinates
            if csys_name != "":
                position = np.dot(transformation, node.coordinates - translation)
            else:
                 position = node.coordinates
            csv_writer.writerow([node.label] + list(position) + list(stress_field[i, :]) + list(strain_field[i, :]) + list(sener_field[i]))

def transform(coordinates, basis1, basis2):
    import numpy as np
