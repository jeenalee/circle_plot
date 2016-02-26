#!/usr/bin/env python

import sys, getopt, cairo, math, numpy

""" 
Draws circle plot for gene expression change between 2 groups of
samples. This code is written for stickleback fish from 2 different
environments, ocean and freshwater.  
"""

fst_file = ""
stat_file = ""
output_dir = ""

myopts, args = getopt.getopt(sys.argv[1:], "hf:s:o:")



def usage():
    print "usage: draw_circle_plot.py:"
    print "  -f <path to fst file>"
    print "  -s <path to statsitcs file>"
    print "  -o <path to output the plot file>"
    sys.exit(1)

    
for o, a in myopts:
    if o == '-h':
        usage()
    if o == '-f':
        fst_file = a
    if o == '-s':
        stat_file = a
    if o == '-o':
        output_dir = a

fst_fh = open(fst_file, "r")
stat_fh = open(stat_file, "r")



# Dictionary 'img' stores all the necessary information for creating
# the size of the output pdf, and defines the center of the circle plot.
img = {}
img['height']     = 1000
img['width']      = 1000
img['center_x']   = img['width']  / 2.0
img['center_y']   = img['height'] / 2.0
img['font_size']  = 16
img['radius_1']   = img['width'] * 0.27
img['radius_2']   = img['width'] * 0.30
img['radius_3']   = img['width'] * 0.31
img['radius_4']   = img['width'] * 0.34 
img['radius_group_label'] = img['width'] * 0.37
img['radius_inner_tick'] = img['width'] * 0.26
img['radius_outer_tick'] = img['width'] * 0.35
img['radius_tick_label'] = img['width'] * 0.245
img['fst_legend_x'] = img['width'] * 0.87
img['stat_legend_x'] = img['width'] * 0.93
img['legend_width'] = 30
img['legend_start_y'] = img['height'] * 0.05
img['legend_height'] = img['legend_width'] * 7


# Global variables to make calling the dictionary easy.
center_x = img['center_x']
center_y = img['center_y']
radius_1 = img['radius_1']
radius_2 = img['radius_2']
radius_3 = img['radius_3']
radius_4 = img['radius_4']
radius_group_label = img['radius_group_label']
radius_inner_tick = img['radius_inner_tick']
radius_outer_tick = img['radius_outer_tick']
radius_tick_label = img['radius_tick_label']
fst_legend_x = img['fst_legend_x'] 
stat_legend_x = img['stat_legend_x']
legend_width = img['legend_width'] 
legend_start_y = img['legend_start_y']
legend_height = img['legend_height']


# The final output file is stored under the path below. The height
# and width of the image can be changed in the img dictionary above.
ps = cairo.PDFSurface(output_dir, img['height'], img['width'])
cr = cairo.Context(ps)
    
cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)


# Tuple to store the number of basepairs in each chromosome.
BP = [('I' , 28185914), ('II' , 23295652), ('III' , 16798506),
      ('IV' , 32632948), ('V' , 12251397), ('VI' , 17083675),
      ('VII' , 27937443), ('VIII' , 19368704), ('IX' , 20249479),
      ('X' , 15657440), ('XI' , 16706052), ('XII' , 18401067),
      ('XIII' , 20083130), ('XIV' , 15246461), ('XV' , 16198764),
      ('XVI' , 18115788), ('XVII' , 14603141), ('XVIII' , 16282716),
      ('XIX' , 20240660), ('XX' , 19732071), ('XXI' ,11717487)]

totalBP = 400788495


# Functions related to parsing data files.

def parse_fst_file(f):
    """
    Parses the fst file without headers and create a dictionary where 
    {'chromosome number' : [(bp, fst), (bp, fst), (bp, fst)]} 
        * fst ranges from 0 to 1.

    """
    fst_dict = {}

    for line in f:
        line = line.split('\t')
        group = line[4]
        
        if "group" in group:
            group = group[5:]
            BP = int(line[5])
            fst = line[24]
            if '-' in fst:
                fst = 0
            fst = float(fst)
            current_entry = (BP, fst)
            
            if group in fst_dict:
                fst_dict[group].append(current_entry)

            else:
                fst_dict[group] = []
                fst_dict[group].append(current_entry)

            fst_dict[group].sort()

    return fst_dict



def parse_stat_file(f):
    """
    Parses the stat file without headers and creates a dictionary:
    {'chromosome number' : [(bp, log(stat)), (bp, log(stat)), 
                            (bp, log(stat))]} 
        * log(stat) ranges from -5 to 5.
    """
    stat_dict = {}

    for line in f:
        line = line.split('\t')
        group = line[1]

        if "group" in group:
            group = group[5:]
            BP = int(line[2])
            stat = math.log(float(line[6]))
        
                
            current_entry = (BP, stat)

                            
            if group in stat_dict:
                stat_dict[group].append(current_entry)

            else:
                stat_dict[group] = []
                stat_dict[group].append(current_entry)

            stat_dict[group].sort()
            
    return stat_dict



# Functions related to finding radian, x and y coordinates, and color.

def get_x_y_coordinates(center_x, center_y, degree, radius):
    """
    Convert a radius and a span of degrees into X, Y coordinates. 
        * Copied from Julian Catchen's slides.
    """
    if degree <= 90:
        theta      = float(degree)
        opp_side   = radius * math.sin(math.radians(theta))
        adj_side   = radius * math.cos(math.radians(theta))
        x          = center_x + adj_side
        y          = center_y + opp_side
    elif degree <= 180:
        theta      = float(degree - 90.0)
        opp_side   = radius * math.sin(math.radians(theta))
        adj_side   = radius * math.cos(math.radians(theta))
        x          = center_x - opp_side
        y          = center_y + adj_side
    elif degree <= 270:
        theta      = float(degree - 180.0)
        opp_side   = radius * math.sin(math.radians(theta))
        adj_side   = radius * math.cos(math.radians(theta))
        x          = center_x - adj_side
        y          = center_y - opp_side
    else:
        theta      = float(degree - 270.0)
        opp_side   = radius * math.sin(math.radians(theta))
        adj_side   = radius * math.cos(math.radians(theta))
        x          = center_x + opp_side
        y          = center_y - adj_side
    return (x, y)


def find_radian(length_to_find_radian, totalBP):
    """
    Given basepair length and total basepair, finds the radian
    appropriate to represent its proportion. Currently, the margin between
    groups are set to 2 degrees.
        * 318.0 = 360.0 - (2.0 * 21 groups)
    """
    radian = 318.0 * (float(length_to_find_radian)/float(totalBP))
    radian = round(radian, 10)
    return radian


def find_color(value_type, value):
    """
    For a given value type (fst or stat), finds the color for the
    value with the shade that represents the value.
    """

    # Transparency is set to 0.8.
    a = 0.8
    
    if value_type == 'fst':
        value = round(float(value), 2)
    
        r = value
        g = value
        b = 0
        
    # Values for representing stat ranges from -5 to 5.     
    elif value_type == 'stat':
        # Multiply by 0.1 so that the values range from -0.5 to 0.5.
        value = round(float(value)*0.1, 2)
        
        # Get the absolute value because fold change is subjective.
        #    It doesn't matter which population expresses gene more
        #    as long as the difference is apparent.
        value = abs(value)
    
        r = 0 + value
        g = 0.5 - value
        b = 0.3

    return (r, g, b, a)




# Functions for drawing arcs and rectangles.

def draw_rectangle(value_type, value, starting_x, starting_y, legend_width, y_increment):
    """ Draws rectangle. To be used for drawing graph legends. """
    # Define colors to use.
    (r, g, b, a) = find_color(value_type, value)
    cr.set_source_rgba(r, g, b, a)

    cr.set_line_width(0.1)

    cr.move_to(starting_x, starting_y)
    cr.line_to(starting_x, starting_y + y_increment)
    cr.line_to(starting_x + legend_width, starting_y + y_increment)
    cr.line_to(starting_x + legend_width, starting_y)
    cr.close_path()
    cr.fill()


def draw_arc(value_type, value, radian_start, radian_end, inner_radius, outer_radius):
    """ Draws arcs given value_type, value, radians and radii. """
    # Define colors to use.
    (r, g, b, a) = find_color(value_type, value)
    cr.set_source_rgba(r, g, b, a)

    cr.set_line_width(0.1)

    # Define starting point on the inner radius.
    (starting_x, starting_y) = get_x_y_coordinates(center_x, center_y, \
                                                   radian_start, inner_radius)
    # Move to starting point.
    cr.move_to(starting_x, starting_y)

    # Draw inner arc.
    cr.arc(center_x, center_y, inner_radius, \
           math.radians(radian_start), math.radians(radian_end))

    # Define second point that is on the outer radius.
    (second_x, second_y) = get_x_y_coordinates(center_x, center_y, radian_end, outer_radius)
    # Draw a straight line to second point.
    cr.line_to(second_x, second_y)
    
    # Draw arc the other direction, and draw the outer arc.
    cr.arc_negative(center_x, center_y, outer_radius, \
                    math.radians(radian_end), math.radians(radian_start))

    # Connect the ends to close the arc.
    cr.close_path()
    cr.fill()

    


# Functions related to annotating the graph.

def write_group_label(group_name, radian_start, radian_end):
    """ Writes group label. """
    cr.set_font_size(15)
    cr.set_source_rgb(0, 0, 0)

    text = group_name
    textents = cr.text_extents(text)
    text_width = textents[2]
    text_height = textents[3]

    radian_group_label = (radian_start + radian_end)/2
    # To adjust the location of the text more center-like, adjust the starting points.
    (group_x, group_y) = get_x_y_coordinates(center_x - (center_x * 0.020), \
                                             center_y + (center_y * 0.015), \
                                             radian_group_label, radius_group_label)

    cr.move_to(group_x, group_y)
    cr.show_text(text)

    
def write_tick_label(tick_name, radian_tick):
    """ Writes tick labels. """
    cr.set_font_size(8)
    cr.set_source_rgb(0.4, 0.5, 0.5)
    
    text = tick_name
    textents = cr.text_extents(text)
    text_width = textents[2]
    text_height = textents[3]

    radian_tick_label = radian_tick - 1
    (tick_x, tick_y) = get_x_y_coordinates(center_x - (center_x * 0.025), \
                                           center_y + (center_y * 0.005), \
                                           radian_tick_label, radius_tick_label)

    cr.move_to(tick_x, tick_y)
    cr.show_text(text)
    
    
def title_legend(value_type):
    """ Writes the title for each legend of colors. """
    legend_annotation_y = legend_start_y - 7
    cr.set_font_size(12)
    cr.set_source_rgb(0, 0, 0)
    
    if value_type == 'fst':
        text = 'Fst'
    if value_type == 'stat':
        text = 'log(FC)'
    
    textents = cr.text_extents(text)
    text_width = textents[2]
    text_height = textents[3]

    if value_type == 'fst':
        cr.move_to(fst_legend_x + 5, legend_annotation_y)
    if value_type == 'stat':
        cr.move_to(stat_legend_x - 3, legend_annotation_y)
        
    cr.show_text(text)

    
def annotate_legend(value_type):
    """ Annotates the scale of colors for the legends. """
    if value_type == 'fst':
        annotation = ['1.0', '0.5', '0.0']
        legend_x = fst_legend_x 
    if value_type == 'stat':
        annotation = ['5.0', '2.5', '0.0']
        legend_x = stat_legend_x 
        
    legend_annotation_x = legend_x - 21
    legend_annotation_y = legend_start_y + 10

    for item in annotation:
        text = item
        textents = cr.text_extents(text)
        text_width = textents[2]
        text_height = textents[3]

        cr.move_to(legend_annotation_x, legend_annotation_y)
        cr.show_text(text)

        legend_annotation_y += (legend_height / 2) - 5

 
    
# Functions related to drawing ticks and legends.

def draw_thicker_line(line_start_x, line_start_y, line_end_x, line_end_y):
    """ Draws thicker line. """
    cr.set_source_rgb(0, 0, 0)
    cr.set_line_width(0.6)
    cr.move_to(line_start_x, line_start_y)
    cr.line_to(line_end_x, line_end_y)
    cr.close_path()
    cr.stroke()

    
def draw_thinner_line(line_start_x, line_start_y, line_end_x, line_end_y):
    """ Draws thinner line. """
    cr.set_source_rgb(0.3, 0.3, 0.3)
    cr.set_line_width(0.3)
    cr.move_to(line_start_x, line_start_y)
    cr.line_to(line_end_x, line_end_y)
    cr.close_path()
    cr.stroke()

    
def make_basepair_ticks(basepair_length, radian_start, tick_length, totalBP):
    """ Draws baspair ticks. """
    individual_tick_radian = find_radian(tick_length, totalBP)

    for i in numpy.arange(0, basepair_length, int(tick_length)):
        tick_number = i / tick_length
        radian_tick = radian_start + (individual_tick_radian * tick_number)
        (tick_start_x, tick_start_y) = get_x_y_coordinates(center_x, center_y, \
                                                           radian_tick, radius_inner_tick)
        (tick_end_x, tick_end_y) = get_x_y_coordinates(center_x, center_y, \
                                                       radian_tick, radius_outer_tick)
        if tick_number % 2 == 0:
            draw_thicker_line(tick_start_x, tick_start_y, tick_end_x, tick_end_y)
            if tick_number !=0 :
                tick_name = str(5 * tick_number) + ' Mb'
                write_tick_label(tick_name, radian_tick)
        else:
            draw_thinner_line(tick_start_x, tick_start_y, tick_end_x, tick_end_y)


def draw_legend(value_type):
    """ Draws the legend of color to indicate the meaning of different shades of colors. """
    value_to_draw = []
    
    if value_type == 'fst':
        start = 0
        end = 1.0
        # Set the increment desired. Too much increments will result
        # in 'stripy' pattern.
        increment = (end - start) / 20
        
        for i in numpy.arange(0, 1.01, increment):
            i = round(i, 2)
            value_to_draw.append(i)

        starting_x = fst_legend_x
        starting_y = legend_start_y
        legend_y_increment = legend_height / len(value_to_draw)
        
        for value in value_to_draw:
            draw_rectangle('fst', value, starting_x, starting_y, \
                           legend_width, legend_y_increment)
            starting_y += legend_y_increment
            
    if value_type == 'stat':
        start = 0.0
        end = 5.0
        # Set the increment desired. Too much increments will result
        # in 'stripy' pattern.
        increment = (end - start) / 20

        # Flip the direction for drawing so that the legend has 0 on
        # the bottom, and 5 on the top.
        for i in numpy.arange(5.01, 0, -increment):
            i = round(i, 2)
            value_to_draw.append(i)

        starting_x = stat_legend_x
        starting_y = legend_start_y
        legend_y_increment = legend_height / len(value_to_draw)
        
        for value in value_to_draw:
            draw_rectangle('stat', value, starting_x, starting_y, \
                           legend_width, legend_y_increment)
            starting_y += legend_y_increment
            
        
        

# Functions related to trimming the data so that it is easier to
# draw. Currently, the basepairs in the data file indicates how far it
# is from the start of the chromosome. These functions (a) find the
# distance between each basepairs of values so that it can find the
# radian based on basepair increment, and (b) bucket basepair increments
# so that small basepair increments can be grouped together to avoid
# drawing arcs that are too small that they overlap. 

def bucket_bp_increment(l_of_t, bucket_threshold):
    """ Buckets basepair increments based on given bucket threshold. """
    bucket_threshold = int(bucket_threshold)
    bp_increment_bucket = []
    interim_storage_bp = 0
    interim_storage_value = 0
    number_of_items_in_interim = 0
    
    for i in range(0, len(l_of_t)):
        bp_increment = l_of_t[i][0]
        value = l_of_t[i][1]

        # First, put the bp_increment and value into the interim storage.
        interim_storage_bp += bp_increment
        interim_storage_value += value
        number_of_items_in_interim += 1

        # If the interim_storage_bp is larger than the bucket, then
        # the interim_storage_bp and the average of values in the
        # interim_storage_value are grouped into a tuple, t.
        # Then, the tuple is appended to the list, bp_increment_bucket.
        if interim_storage_bp >= bucket_threshold:

            avg_value = interim_storage_value / number_of_items_in_interim
            t = (interim_storage_bp, avg_value)
            bp_increment_bucket.append(t)

            # The interims are reset to 0.
            interim_storage_bp = 0
            interim_storage_value = 0
            number_of_items_in_interim = 0

    # After going through every item in the input list of tuples,
    # l_of_t, if there's any item left in the interim storage, then
    # add the bp_increment and average value to the tuple, t, and then
    # to the list, bp_increment_bucket.
    if number_of_items_in_interim > 0:
        avg_value = interim_storage_value / number_of_items_in_interim
        t = (interim_storage_bp, avg_value)
        bp_increment_bucket.append(t)

    return bp_increment_bucket
            

def find_bp_increment(l_of_t, group_length):
    """
    Finds distance between a pair of basepair locations. The last
    basepair location is stretched until the end of the chromosome.
    """
    l_bp_increment = []
    for i in range(0, len(l_of_t)):

        # The l_of_t is in the following format:
        # l = [(bp_from_start_of_chromosome, some_fst_or_stat_value), (bp2, value2), ...]
        value = float(l_of_t[i][1])

        # For the first tuple in list, bp_increment would remain the same.
        if i == 0:
            bp_increment = int(l_of_t[i][0])

        # For the last basepair location, the second to last bp
        # location is subtracted from the length of the chromosome. By
        # doing so, the value associated with the last basepair
        # location is applied to the end of the chromosome.
        elif i == len(l_of_t):
            bp_increment = int(group_length) - int(l_of_t[i-2][0])

        # Otherwise, the bp_increment is calculated by subtracting the
        # previous basepair location from the current basepair
        # location.
        else:
            bp_increment = int(l_of_t[i][0]) - int(l_of_t[i-1][0])

        # The bp_increment and value are added to a tuple, t.
        t = (bp_increment, value)
        # The tuple t is added to the list, l_bp_increment.
        l_bp_increment.append(t)

    # l_bp_increment is bucketed through the function
    # bucket_bp_increment. You can change the bucket threshold by
    # specifying it below. Too large bucket size results in "chunky"
    # looking circle plot. Too small bucket size results in drawing
    # arcs that are too small that they overlap, creating a "stripy"
    # looking circle plot.
    
    l_bp_increment = bucket_bp_increment(l_bp_increment, 300000)
    return l_bp_increment




# Last set of functions to call previous functions and draw circle
# plot, and add legends.

def draw_circle_plot(l_of_t, totalBP, fst_f, stat_f):
    """ Draws circle plot. """

    # The margin between groups is 2 degrees.
    radian_start = 2.0
    radian_end = 0
    
    fst_dict = parse_fst_file(fst_f)
    stat_dict = parse_stat_file(stat_f)

    for (group_name, group_length) in l_of_t:

        make_basepair_ticks(group_length, radian_start, 5000000, totalBP)

        radian_start_fst = radian_start
        radian_start_stat = radian_start
        radian_end_fst = radian_end
        radian_end_stat = radian_end
        
        fst_l = fst_dict[group_name]
        stat_l = stat_dict[group_name]
        
        fst_l_by_increment = find_bp_increment(fst_l, group_length)
        stat_l_by_increment = find_bp_increment(stat_l, group_length)

        for (bp, value) in fst_l_by_increment:
            angle_for_bp_fst = find_radian(bp, totalBP)
            radian_end_fst = angle_for_bp_fst + radian_start_fst
            draw_arc('fst', value, radian_start_fst, radian_end_fst, radius_3, radius_4)
            radian_start_fst = radian_end_fst
        
        for (bp, value) in stat_l_by_increment:
            
            angle_for_bp_stat = find_radian(bp, totalBP)
            radian_end_stat = angle_for_bp_stat + radian_start_stat
            draw_arc('stat', value, radian_start_stat, radian_end_stat, radius_1, radius_2)
            radian_start_stat = radian_end_stat


        radian_end = radian_end_fst
        
        write_group_label(group_name, radian_start, radian_end)
        radian_start = (radian_end + 2.0)
        

def add_legends(legend1, legend2):
    """ Adds legends. """
    draw_legend(legend1)
    title_legend(legend1)
    annotate_legend(legend1)

    draw_legend(legend2)
    title_legend(legend2)
    annotate_legend(legend2)
    



# Now calling functions!

draw_circle_plot(BP, totalBP, fst_fh, stat_fh)
add_legends('fst', 'stat')


cr.show_page()
