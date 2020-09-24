#author: siayouyang
#two-component multivalent interactions simulation
#refer to P., Li, 2012, nature.

from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import random
import math

#OptionParser
parser = OptionParser()
parser.add_option("--m1", dest="molecule1", default=20, type=int, help="number of molecule 1 (default=20)")
parser.add_option("--m2", dest="molecule2", default=20, type=int, help="number of molecule 2 (default=20)")
parser.add_option("-f", "--frames", dest="updates", default=10000, type=int, help="number of frames to be updated (default=10000)")
parser.add_option("-v", "--valencies", dest="domains", default=4, type=int, help="number of domains in a molecule (valency) (default=4)")
parser.add_option("-d", "--dimensions", dest="dimensions", default=10, type=int, help="define the volume of a whole cube (default=10)")
parser.add_option("-o", "--output", dest="output_name", default="phase_valency.csv", help="output name, recommend csv format (default=phase_valency.csv)")
(options, args) = parser.parse_args()

fig = plt.figure()
ax = p3.Axes3D(fig)
molecule1 = options.molecule1    #set number of molecule 1
molecule2 = options.molecule2   #set number of molecule 2
updates = options.updates  #set updates(frames)
domains = options.domains    # set the number of domains in each element(molecule)
dimensions = options.dimensions  #eg. dimensions = 10 means 10*10*10 units cube
output_name = options.output_name

'''''''''
index information and location:
[x,y,z,size,red,blue,xyz,minusplus, accred, accblue, free_domain, total_molecule, status, label,domain...]
x = data[iteration][j, 0:1]
y = data[iteration][j, 1:2]
z = data[iteration][j, 2:3]
size = data[iteration][j, 3:4]
red = data[iteration][j, 4:5]
blue = data[iteration][j, 5:6]
xyz = data[iteration][j, 6:7]
minusplus = data[iteration][j, 7:8]
accred = data[iteration][j, 8:9]
accblue = data[iteration][j, 9:10]
free_domain = data[iteration][j, 10:11]
total_molecule = data[iteration][j, 11:12] (total molecules 1 or 2 sharing identical coordinate/cube)
status = data[iteration][j, 12:13] (1=less, 2=more, 3=equal, 4=no partner)
label = data[iteration][j, 13:14]
domain = data[iteration][j, 14:]
'''''''''
# write and save as csv
csv = open(output_name, 'w')
csv.write("t" + ",count" + ",max_complex\n")

# define number of indices according to number of domains
dims = (int(14 + domains) , 1)
initial_position_list_a = [[0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0]]
for i in range(domains):
    initial_position_list_a.append([0])

initial_position_list_b = [[dimensions], [dimensions], [dimensions], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1]]
for i in range(domains):
    initial_position_list_b.append([1])

# Random initial positions.
start_positions = np.array(list(
    map(np.random.randint, initial_position_list_a, initial_position_list_b, [molecule1+molecule2] * dims[0]))).T

data = [start_positions]
# define two color/two different molecules
for init_red in range(molecule1):
    data[0][init_red, 4:5] = data[0][init_red, 4:5] + 1

for init_blue in range(molecule2):
    data[0][init_blue+molecule1, 5:6] = data[0][init_blue+molecule1, 5:6] + 1

# iteration starts from here!
particles_list = []
max_complex_list = []
for iteration in range(updates+1):
    # reset indices(free_domain, total_elements, status, label)
    for j in range(data[0].shape[0]):
        data[iteration][j, 10:14] = [0, 0, 0, 0]

    # dissociate
    for j in range(data[0].shape[0]):
        for d in range(len(data[iteration][j, 14:])):
            if int(data[iteration][j, (14+d):(15+d)]) != 0:
                yes_or_no = random.choices([1, 2], weights = [1, 49], k = 1)     # weights : dissociation probability/2 = (0.1*(1-0.6))/2 = 0.04/2 = 0.02
                if yes_or_no == [1]:   #yes, they dissociate
                    n = int(data[iteration][j, (14 + d):(15 + d)])
                    for d2 in range(len(data[iteration][(n-1), 14:])):
                        if int(data[iteration][(n-1), (14 + d2):(15 + d2)]) == (j+1):
                            data[iteration][(n-1), (14 + d2):(15 + d2)] = [0]    # remove its domain from its binding partner
                            break
                        else:
                            pass
                    data[iteration][j, (14 + d):(15 + d)] = [0]                  # remove its partner's domain from itself
                elif yes_or_no == [2]:   #No, they retain associated
                    pass
            elif int(data[iteration][j, (14 + d):(15 + d)]) == 0:
                pass

    # update free_domain
    for j in range(data[0].shape[0]):
        temp_zero = []
        for d in range(len(data[iteration][j, 14:])):
            if int(data[iteration][j, (14+d):(15+d)]) == 0:
                zero = int(data[iteration][j, (14 + d):(15 + d)])
                temp_zero.append(zero)
            else:
                pass
        data[iteration][j, 10:11] = len(temp_zero)

    # update total_element
    for j in range(data[0].shape[0]):
        for k in range(data[0].shape[0]):
            if data[iteration][j, 0:3].tolist() == data[iteration][k, 0:3].tolist() and data[iteration][j, 4:6].tolist() == data[iteration][k, 4:6].tolist():
                data[iteration][k, 11:12] = data[iteration][k, 11:12] + 1
            else:
                pass

    # update status(less/more/equal/no partner)
    for j in range(data[0].shape[0]):
        for k in range(data[0].shape[0]):
            if data[iteration][j, 0:3].tolist() == data[iteration][k, 0:3].tolist() and data[iteration][j, 4:6].tolist() != data[iteration][k, 4:6].tolist():
                if data[iteration][j, 11:12] < data[iteration][k, 11:12]:
                    data[iteration][j, 12:13] = 1  # 1=less
                    break
                elif data[iteration][j, 11:12] > data[iteration][k, 11:12]:
                    data[iteration][j, 12:13] = 2  # 2=more
                    break
                elif data[iteration][j, 11:12] == data[iteration][k, 11:12]:
                    data[iteration][j, 12:13] = 3  # 3=equal
                    break
            else:
                pass  # 0=no partner(default)


    def association():
        if data[iteration][j, 10:11] != [0]:  # there are free domains left
            for d in range(len(data[iteration][j, 14:])):
                if int(data[iteration][j, (14 + d):(15 + d)]) == 0:  # location of a free domain is found
                    temp_counterparts = []
                    for k in range(data[0].shape[0]):  # searching for counterparts (with different color)
                        if data[iteration][j, 0:3].tolist() == data[iteration][k, 0:3].tolist() and data[iteration][j, 4:6].tolist() != data[iteration][k, 4:6].tolist():
                            for d2 in range(len(data[iteration][k, 14:])):  # searching for counterparts' free domain
                                if int(data[iteration][k, (14 + d2):(15 + d2)]) == 0:
                                    temp_counterparts.append(k)  # append to a list for later lottery
                                else:
                                    pass
                    random.shuffle(temp_counterparts)  # colide to each domain randomly without repetition until all tried or bound
                    for l in range(len(temp_counterparts)):
                        total_in_cube = int(data[iteration][j, 11:12]) + int(data[iteration][temp_counterparts[l], 11:12])  # reflect the crowdedness(total number of domain in a cube)
                        association_prob = round(((0.6 * 0.9 * 20) / (20 + total_in_cube * domains)) * 1000)  # ensure thet there is 4 significant digits
                        yes_or_no = random.choices([1, 2], weights=[association_prob, (1000 - association_prob)], k=1)  # association probability
                        if yes_or_no == [1]:  # yes, they associate
                            if data[iteration][j, (14 + d):(15 + d)] == 0:
                                data[iteration][j, (14 + d):(15 + d)] = [temp_counterparts[l] + 1]  # record their connections
                                data[iteration][j, 10:11] = data[iteration][j, 10:11] - 1  # update the number of free domain left
                                for d3 in range(len(data[iteration][temp_counterparts[l], 14:])):  # searching for counterparts' free domain
                                    if int(data[iteration][temp_counterparts[l], (14 + d3):(15 + d3)]) == 0:
                                        data[iteration][temp_counterparts[l], (14 + d3):(15 + d3)] = [j + 1]  # record their connections
                                        data[iteration][temp_counterparts[l], 10:11] = data[iteration][temp_counterparts[l], 10:11] - 1  # update the number of free domain left
                                        break
                                    else:
                                        pass
                            else:
                                break
                        elif yes_or_no == [2]:
                            pass

    # associate
    for j in range(data[0].shape[0]):
        if data[iteration][j, 12:13] == 1:   # using the fewer elements' domains to cycles through theirs counterparts
            association()
        elif data[iteration][j, 12:13] == 3 and data[iteration][j, 4:6].tolist() == [1, 0]:   # for equal elements, red domains will be used to cycles through theirs counterparts
            association()

    #updata label (elements/molecules with same label will diffuse together)
    for j in range(data[0].shape[0]):
        for k in range(data[0].shape[0]):
            if data[iteration][j, 0:3].tolist() == data[iteration][k, 0:3].tolist() and j != k:
                k_elements = (data[iteration][k, 14:]).tolist()
                for n in range(len(k_elements)):
                    if 0 in k_elements:
                        k_elements.remove(0)
                    else:
                        pass
                if (j+1) in k_elements:              # if there is a connection between them, update their label
                    if int(data[iteration][j, 13:14]) == 0 and int(data[iteration][k, 13:14]) == 0:
                        data[iteration][j, 13:14] = [j+1]
                        data[iteration][k, 13:14] = [j+1]
                    elif int(data[iteration][j, 13:14]) != 0 and int(data[iteration][k, 13:14]) == 0:
                        data[iteration][k, 13:14] = data[iteration][j, 13:14]
                    elif int(data[iteration][j, 13:14]) == 0 and int(data[iteration][k, 13:14]) != 0:
                        data[iteration][j, 13:14] = data[iteration][k, 13:14]
                    elif int(data[iteration][j, 13:14]) != 0 and int(data[iteration][k, 13:14]) != 0:
                        if int(data[iteration][j, 13:14]) > int(data[iteration][k, 13:14]):
                            for l in range(data[0].shape[0]):
                                if int(data[iteration][l, 13:14]) == int(data[iteration][j, 13:14]) and l != j:
                                    data[iteration][l, 13:14] = data[iteration][k, 13:14]
                                else:
                                    pass
                            data[iteration][j, 13:14] = data[iteration][k, 13:14]
                        elif int(data[iteration][j, 13:14]) < int(data[iteration][k, 13:14]):
                            for l in range(data[0].shape[0]):
                                if int(data[iteration][l, 13:14]) == int(data[iteration][k, 13:14]) and l != k:
                                    data[iteration][l, 13:14] = data[iteration][j, 13:14]
                                else:
                                    pass
                            data[iteration][k, 13:14] = data[iteration][j, 13:14]
                        elif int(data[iteration][j, 13:14]) == int(data[iteration][k, 13:14]):
                            pass

    # if they have the same label they diffuse together
    for j in range(data[0].shape[0]):
        if data[iteration][j, 13:14] == 0:                 # molecules with no partner will 100% diffuse
            temp_xyz = np.random.choice([1, 2, 3])
            temp_minusplus = np.random.choice([-1, 1])
            data[iteration][j, 6:7] = temp_xyz
            data[iteration][j, 7:8] = temp_minusplus
        elif data[iteration][j, 13:14] != 0:
            temp_molecules = []
            for k in range(data[0].shape[0]):               # molecules with same coordinate and label(not 0)
                if data[iteration][j, 0:3].tolist() == data[iteration][k, 0:3].tolist() and data[iteration][j, 13:14] == data[iteration][k, 13:14]:
                    temp_molecules.append(k)
                else:
                    pass
            total_molecules = len(temp_molecules)
            diffuse_probability = round((1/math.sqrt(total_molecules))*1000)
            remove_or_not = random.choices([1, 2], weights=[diffuse_probability, (1000-(diffuse_probability))], k=1)
            if remove_or_not == [1]:
                temp_xyz = np.random.choice([1, 2, 3])
                temp_minusplus = np.random.choice([-1, 1])
            elif remove_or_not == [2]:
                temp_xyz = 0
                temp_minusplus = 0
            for k in range(data[0].shape[0]):
                if data[iteration][j, 0:3].tolist() == data[iteration][k, 0:3].tolist() and data[iteration][j, 13:14] == data[iteration][k, 13:14]:
                    data[iteration][j, 6:7] = temp_xyz
                    data[iteration][k, 6:7] = data[iteration][j, 6:7]
                    data[iteration][j, 7:8] = temp_minusplus
                    data[iteration][k, 7:8] = data[iteration][j, 7:8]
                else:
                    pass

    previous_data = data[iteration]
    # Movement
    for j in range(data[0].shape[0]):
        xyz = int(data[iteration][j, 6:7])
        minusplus = int(data[iteration][j, 7:8])

        # confine the volume
        if xyz == 0 and minusplus == 0:
            pass
        else:
            box = np.arange(0,dimensions,1)
            if int(previous_data[j, (xyz-1):xyz] + minusplus) in box:
                previous_data[j, (xyz-1):xyz] = previous_data[j, (xyz-1):xyz] + minusplus
            else:
                pass

    # count particles left
    data_array = []
    for j in range(data[0].shape[0]):
        list = (data[iteration][j, 0:3]).tolist()
        if list not in data_array:
            data_array.append(list)
    particles = len(data_array)
    particles_list.append(particles)

    # number of molecules in a single complex
    data_array_2 = []
    for j in range(data[0].shape[0]):
        temp_count_molecules = []
        for k in range(data[0].shape[0]):
            if data[iteration][j, 0:3].tolist() == data[iteration][k, 0:3].tolist() and data[iteration][j, 13:14] == data[iteration][k, 13:14]:
                temp_count_molecules.append(k)
        count_molecules = len(temp_count_molecules)
        data_array_2.append(count_molecules)
    max_complex = max(data_array_2)
    max_complex_list.append(max_complex)

    # append data of current iteration
    data.append(np.array(previous_data))

# let the droplets in the same coordinate become larger
for i in range(updates + 1):
    for j in range(data[0].shape[0]):
        for k in range(data[0].shape[0]):
            if data[i][j, 0:3].tolist() == data[i][k, 0:3].tolist():
                data[i][j, 3:4] = data[i][j, 3:4] + 1
            else:
                pass

# mixing colors
for i in range(updates + 1):
    for j in range(data[0].shape[0]):
        for k in range(data[0].shape[0]):
            if data[i][j, 0:3].tolist() == data[i][k, 0:3].tolist():
                data[i][j, 8:9] = data[i][j, 8:9] + int(data[i][k, 4:5])
                data[i][j, 9:10] = data[i][j, 9:10] + int(data[i][k, 5:6])
            else:
                pass

# Provide starting angle for the view.
ax.view_init(20, 45)


def scatters(j):
    #clear the previous scatters
    ax.clear()
    # Setting the axes properties
    ax.set_xlim3d([0, dimensions])
    ax.set_xlabel('X')

    ax.set_ylim3d([0, dimensions])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0, dimensions])
    ax.set_zlabel('Z')

    ax.set_title('3D Animated Phase Saperation')

    #update scatter plot
    for k in range(data[0].shape[0]):
        accred = float(data[j][k, 8:9])
        accblue = float(data[j][k, 9:10])
        normred = float(accred / (accred + accblue))
        normblue = float(accblue / (accred + accblue))
        ax.scatter(data[j][k, 0:1], data[j][k, 1:2], data[j][k, 2:3], color=(normred, 0, normblue), s=30 * (data[j][k, 3:4]))
    plt.suptitle("t = %d" % j + "    particles = %d" % particles_list[j] + "    the largest has %d molecules" % max_complex_list[j])
    csv.write("%d"%j + ",%d"% particles_list[j] + ",%d\n"% max_complex_list[j])


def animate(j):
    scatters(j+1)
    return scatters


def init():
    scatters(0)
    return scatters

ani = animation.FuncAnimation(fig=fig, func=animate, frames=updates, init_func=init, interval=10, blit=False, repeat=False)

plt.show()
csv.close()