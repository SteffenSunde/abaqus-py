def field_to_csv(odb_handle, elset_handle, field_name, step_number=-1, frame_number=-1):
    from abaqusConstants import CENTROID
    import csv

    frame_handle = odb_handle.steps.items()[step_number][1].frames[frame_number]
    field = frame_handle.fieldOutputs[field_name].getSubset(region=elset_handle, position=CENTROID)
    file_name = 'CENTROID-{}-{}.csv'.format(elset_handle.name, field.name)
    with open(file_name, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')
        csv_writer.writerow(('Step time:{}'.format(frame_handle.frameValue),))
        csv_writer.writerow(('Element_id',) + field.componentLabels)
        for point in field.values:
            csv_writer.writerow((point.elementLabel,) + tuple(point.data))


def field_history_to_csv(odb_handle, elset_handle, field_name, steps=[1], frequency=1):
    # Create folder?
    for step in steps:
        step = odb_handle.steps.items()[step][1]
        if frequency < len(step.frames) and frequency > 1:
            frames = distributed_selection(step.frames, frequency)
            #frames = [step.frames[-1]]
            for frame in frames:
                field_to_csv(odb_handle, elset_handle, field_name)


def distributed_selection(lst, num):
    """
    Selects evenly elements from the given list with the given number of elements.
    Ensures that the last point are included. The first point is not necessary due
    to being part of the previous time step.

    For two points, the mid and last point should be selected...
    
    Not thuroughly tested

    """

    if num > len(lst):
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

