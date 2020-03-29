from abaqusConstants import (CENTROID, SCALAR, ELEMENT_NODAL, NODAL, DEFORMED, 
    TRUE_DISTANCE, PATH_POINTS, OFF)
import csv
import numpy as np
import operator

def field_to_csv(odb_handle, elset_handle, field_name, step_number=-1, frame_number=-1):
    """
    Outputs the given field at given frame number in odb file.

    TODO:
    - Make similar for nodal quantities
    - Test
    """
    
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


def distributed_selection(lst, num, include_end=True):
    """
    Selects evenly elements from the given list with the given number of elements.
    Ensures that the last point are included. The first point is not necessary due
    to being part of the previous time step.

    For two points, the mid and last point should be selected...
    
    TODO:
    - Test
    - Abaqus sometimes have a really small last increment to meet time step.
      in these cases, maybe skipping the second last element is better than the first?
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

    # Depends on the field being ordered with ascending node label? No
    set_handle = odb_handle.rootAssembly.nodeSets[set_name.upper()]
    if csys_name == "":
        field = frame_handle.fieldOutputs[field_name].getSubset(region=set_handle, position=ELEMENT_NODAL)
    else:
        deformation_field = frame_handle.fieldOutputs['U']
        csys_handle = odb_handle.rootAssembly.datumCsyses[csys_name.upper()]
        field = frame_handle.fieldOutputs[field_name].getSubset(region=set_handle, position=ELEMENT_NODAL).getTransformedField(datumCsys=csys_handle, deformationField=deformation_field)

    num_components = max(len(field.componentLabels), 1)
    result_dict = {}
    connected_dict = dict([(node.nodeLabel, 0) for node in field.values])

    # This doesnt work because for node set 
    # SLS 28.03: Seems to give half the correct value..
    for element_node in field.values: 
        node_label = element_node.nodeLabel  
        # result[i, :] += element_node.data
        # connected_elements[i] += 1
        if node_label in result_dict:
            result_dict[node_label] += element_node.data
        else:
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
    # TODO: Doesnt work for un-even spaced nodes (mid-side nodes)!
    result = np.ndarray((len(result_dict), num_components))
    for i, node in enumerate(set_handle.nodes[0]):
        node_label = node.label
        result[i, :] = np.asarray(result_dict[node_label]) / connected_dict[node_label]
        # result[i, -1] = node_label  # TODO

    return result, labels


def create_rotating_csys(odb_handle, origin, x, y, name):
    node1 = odb_handle.rootAssembly.nodeSets[origin].nodes[0][0]
    node2 = odb_handle.rootAssembly.nodeSets[x].nodes[0][0]
    node3 = odb_handle.rootAssembly.nodeSets[y].nodes[0][0]
    odb_handle.rootAssembly.DatumCsysByThreeNodes(name=name, 
        coordSysType=CARTESIAN, origin=node1, point1=node2, point2=node3)

    if odb_handle.isReadOnly():
        print("Odb is read-only, so csys {} is only stored for the current session".format(name))
    else:
        odb_handle.save()
    return

def contact_data(field_name, odb_handle, set_name, frame_handle, csys_name):
    """
    Extract contact information from odb file..

    """
    # Depends on the field being ordered with ascending node label? No
    set_handle = odb_handle.rootAssembly.nodeSets[set_name.upper()]
    if csys_name == "":
        field = frame_handle.fieldOutputs[field_name].getSubset(region=set_handle, position=NODAL)
    else:
        csys_handle = odb_handle.rootAssembly.datumCsyses[csys_name.upper()]
        field = frame_handle.fieldOutputs[field_name].getSubset(region=set_handle, position=NODAL).getTransformedField(datumCsys=csys_handle)

    num_components = max(len(field.componentLabels), 1)
    # result_dict = {}
    # connected_dict = dict([(node.nodeLabel, 0) for node in field.values])

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


def nodal_data_to_csv_section(odb_handle, set_name, step_number=-1, frame_number=-1, csys_name=''):
    """
    TODO: Check SENER output!
    TODO: Fix for arbitrary node set!

    """

    if csys_name != "":
        csys_handle = odb_handle.rootAssembly.datumCsyses[csys_name.upper()]
        transformation =    np.outer(csys_handle.xAxis, [1, 0, 0]) + \
                            np.outer(csys_handle.yAxis, [0, 1, 0]) + \
                            np.outer(csys_handle.zAxis, [0, 0, 1])
        transformation.transpose()
        #translation = csys_handle.origin  # What if csys is attached to node? origin will be [0, 0, 0]
        #translation = odb_handle.rootAssembly.nodeSets[set_name.upper()].instances[0].nodes[1349].coordinates
        translation = odb_handle.rootAssembly.nodeSets['NODEORIGIN'].nodes[0][0].coordinates
        print("CSYS with translation {} chosen.".format(translation))

    # nodal_data_section(field_name, odb_handle, set_name, frame_handle):
    #instance_handle = odb_handle.rootAssembly.instances[part_name.upper()]
    set_handle = odb_handle.rootAssembly.nodeSets[set_name.upper()]
    step_handle =  odb_handle.steps.items()[step_number][1]
    frame_handle = step_handle.frames[frame_number]
    job_name = odb_handle.name.split('/')[-1][:-4]
    file_name = '{}_NODAL-{}-{}-{}.csv'.format(job_name, set_name, step_number, frame_number)
    stress_field, stress_labels = nodal_data_section('S', odb_handle, set_name, frame_handle, csys_name)
    strain_field, strain_labels = nodal_data_section('LE', odb_handle, set_name, frame_handle, csys_name)
    sener_field, sener_label = nodal_data_section('SENER', odb_handle, set_name, frame_handle)
    total_time = step_handle.totalTime + frame_handle.frameValue

    # TODO: Add for CPRESS etc. 
    # For coordinates:
    # Find translation?
    # Rotate?

    print("Exporting {}".format(file_name))
    with open(file_name, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', lineterminator = '\n')
        csv_writer.writerow(['Step time:{}'.format(total_time)])
        csv_writer.writerow(['Node_label', 'x', 'y', 'z'] + stress_labels + strain_labels + sener_label)
        #csv_writer.writerow(['Node_label', 'x', 'y', 'z'] + stress_labels + strain_labels)
        for i, node in enumerate(set_handle.nodes[0]):
            node_label = node.label
            #position = instance_handle.nodes[node_label - 1].coordinates
            if csys_name != "":
                position = np.dot(transformation, node.coordinates - translation)
            else:
                 position = node.coordinates
            csv_writer.writerow([node.label] + list(position) + list(stress_field[i, :]) + list(strain_field[i, :]) + list(sener_field[i]))
            # csv_writer.writerow([node.label] + list(position) + list(stress_field[i, :]) + list(strain_field[i, :]))
    # return i


def extract_history(odb_handle, set_name, increment_list, csys_name=''):
    """
    Extracts the data for the list of increments given. Extracts for all available increments
    if step contains less than the specified number of increments.

    TODO: When time series is extracted, first of one time step and the last of the next time 
    step have the same time value. Therefore, add one to the k chosen elements and remove the
    last?

    """
    num_nodes = len(odb_handle.rootAssembly.nodeSets[set_name.upper()].nodes[0])
    total_steps = 0
    for i, inc in enumerate(increment_list):
        increments = odb_handle.steps[odb_handle.steps.items()[i][1].name].frames
        max_increments = len(increments)
        if inc > max_increments:
            print("Warning: Maximum number of frames from step {} is {}. ({} was requested).".format(i, max_increments, inc))
            increments = [frame.frameId for frame in increments]
        elif inc == 0:
            increments = []
        else:
            hash_map = dict([(frame.frameValue, frame.frameId) for frame in increments])
            increments = choose_subset([frame.frameValue for frame in increments], inc+1)
            increments = [hash_map[k] for k in increments[1:]]

        for j in increments:
            nodal_data_to_csv_section(odb_handle, set_name, i, j, csys_name)
            total_steps += 1

    print("A total of {} time steps was extracted successfully for {} nodes.".format(total_steps, num_nodes))
    return


def extract_history_adjusted(odb_handle, set_name, increment_list, csys_name=''): # TODO: Find better name
    # Not working? See function extract_history().
    for i, inc in enumerate(increment_list):
        step_handle = odb_handle.steps[odb_handle.steps.items()[i][1].name]
        increments = step_handle.frames
        max_increments = len(increments)
        if inc > max_increments:
            increments = [frame.frameId for frame in increments]
        elif inc == 0:
            increments = []
        else:
            increments = [frame.frameValue for frame in increments]
            hash_map = dict(zip(range(len(increments)), increments))
            even_spread = choose_subset(increments, inc + 1)  # +1 because method includes both ends.

            increments = [hash_map[x] for x in even_spread]
            print(increments)

        for j in increments:
            nodal_data_to_csv_section(odb_handle, set_name, csys_name=csys_name, step_number=i, frame_number=j, )

    return


def transform(coordinates, basis1, basis2):
    return


def extract_contact(odb_handle, set_name, step_name, csys_name, frame_no=-1, meta=''):
    """
    Extracts the following fields from a given contact node set:
        x, y, cpress, cshear1, cslip1, s11 (rotated).
    
    The nodes are unordered in the node set! Therefore, data is sorted according
    to increasing x-values before exported. This means that the node set is assumed
    to be oriented from left to right.
    """

    # Pointer to abaqus objects
    step_handle = odb_handle.steps[step_name]
    step_no = step_handle.number
    frame_handle = step_handle.frames[frame_no]
    csys_handle = odb_handle.rootAssembly.datumCsyses[csys_name]
    set_handle = odb_handle.rootAssembly.nodeSets[set_name.upper()]
    num_nodes = len(set_handle.nodes[0])
    instance_handle = odb_handle.rootAssembly.instances[set_handle.instanceNames[0]]

    # Allocate space
    data = np.ndarray((num_nodes, 6))

    # Retrieve fields from odb
    cpress = frame_handle.fieldOutputs['CPRESS'].getSubset(region=set_handle, position=NODAL)
    cshear = frame_handle.fieldOutputs['CSHEAR1'].getSubset(region=set_handle, position=NODAL)
    cslip = frame_handle.fieldOutputs['CSLIP1'].getSubset(region=set_handle, position=NODAL)
    deformation_field = frame_handle.fieldOutputs['U']
    s11 = frame_handle.fieldOutputs['S'].getSubset(region=set_handle, position=ELEMENT_NODAL).getTransformedField(datumCsys=csys_handle, deformationField=deformation_field)

    # Interpolate nodal values for stress from element_nodal values.
    node_map = {}
    visits = {}
    for value in s11.values:
        node_label = value.nodeLabel
        if node_label in node_map:
            node_map[node_label] += value.data[0]
            visits[node_label] += 1
        else:
            node_map[node_label] = value.data[0]
            visits[node_label] = 1

    s11 = {key: node_map[key]/visits[key] for (key, _) in node_map.iteritems()}

    # Fill in numpy ndarray
    for i in range(num_nodes):
        node_label = cpress.values[i].nodeLabel
        coords = instance_handle.nodes[node_label - 1].coordinates
        data[i, 0] = coords[0]
        data[i, 1] = coords[1]
        data[i, 2] = cpress.values[i].data
        data[i, 3] = cshear.values[i].data
        data[i, 4] = cslip.values[i].data
        data[i, 5] = s11[node_label]
    
    # Sort entries in data according to x (should then sort by y?)
    data = data[data[:,0].argsort()]

    # Export data
    # header argument in savetxt() is not available in abaqus version of numpy btw
    file_name = "contact_parameters-{}-{}-{}.csv".format(step_no, frame_no, meta)
    with open(file_name, 'w') as f:
        f.write("x,y,cpress,cshear1,cslip1,s11\n")
    with open(file_name, 'a') as f:
        np.savetxt(f, data, delimiter=',')
        print("Data successfully written to {}.".format(file_name))


def extract_contact_path(odb_handle, session_handle, path_name, steps = [], meta=''):
    """
    Remember when reading csv files using pandas 
    library, use delim_whitespace=True
    TODO: Add extraction of quantities for Ruiz parameters
    """

    if not steps:
        steps = [odb_handle.steps[step_name].number-1 for step_name,_ in odb_handle.steps.items()]

    pth = session_handle.paths[path_name]
    session_handle.viewports['Viewport: 1'].setValues(displayedObject=odb_handle)

    for step in steps:
        session_handle.viewports['Viewport: 1'].odbDisplay.setFrame(step=step, frame=-1)

        session_handle.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
            variableLabel='CPRESS', outputPosition=ELEMENT_NODAL)
        session_handle.XYDataFromPath(name='cpress', path=pth, includeIntersections=False, 
            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
            projectionTolerance=0, shape=DEFORMED, labelType=TRUE_DISTANCE)

        session_handle.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
            variableLabel='CSLIP1', outputPosition=ELEMENT_NODAL)
        session_handle.XYDataFromPath(name='cshear', path=pth, includeIntersections=False, 
            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
            projectionTolerance=0, shape=DEFORMED, labelType=TRUE_DISTANCE)

        session_handle.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
            variableLabel='CSLIP1', outputPosition=ELEMENT_NODAL)
        session_handle.XYDataFromPath(name='cslip', path=pth, includeIntersections=False, 
            projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, 
            projectionTolerance=0, shape=DEFORMED, labelType=TRUE_DISTANCE)

        x0 = session_handle.xyDataObjects['cpress']
        x1 = session_handle.xyDataObjects['cshear']
        x2 = session_handle.xyDataObjects['cslip']
        
        session_handle.writeXYReport(fileName='contact-{}-{}.csv'.format(step, meta), appendMode=OFF, xyData=(x0, x1, x2))
        
##########################################################

# Abaqus doesnt like imports to have dependencies. 

def choose_subset(numbers, k, sorted = True):
    """Chooses k elements from the input list of positive numbers, maximising scatter."""

    # Edge cases
    n = len(numbers)
    if k >= n:
        return numbers
    elif k < 0:
        return []

    # Sort input array if it isn't already
    if not sorted:
        numbers.sort()

    # Allocate the dynamic programming (dp) table and parent pointers (pp)
    dp = -1*np.ones((k, n))
    pp = np.zeros((k, n))

    # Edge case, because choosing two points from a sorted list is always the extreme
    for j in range(1, n):
        dp[1, j] = numbers[j] - numbers[0]
        pp[1, j] = 0
    
    # Fill in rest of the DP table
    for i in range(2, k):
        for j  in range(i, n):
            candidates = []
            for l in range(0, j):
                candidates.append(min(dp[i-1, l], numbers[j] - numbers[l]))
            pp[i, j], dp[i, j] = max(enumerate(candidates), key=operator.itemgetter(1))

    # Backtrace the DP table to choose the k optimal points
    result = [numbers[-1]]
    j = n-1
    for i in range(k-1, 0, -1):
        result.append(numbers[int(pp[i, j])])
        j = pp[i, j]
    result.reverse()
    
    return result
