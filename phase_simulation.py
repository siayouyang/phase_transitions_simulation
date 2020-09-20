#author: siayouyang
#two-component multivalent interactions simulation
#refer to P., Li, 2012, nature.

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

fig = plt.figure()
ax = p3.Axes3D(fig)
elements = 50    #set elements
updates = 100   #set updates

# [x,y,z,size,red,blue,xyz,minusplus, accred, accblue]
dims = (10, 1)

# Random initial positions.
start_positions = np.array(list(
    map(np.random.randint, [[0], [0], [0], [0], [0], [0], [0], [0], [0], [0]], [[10], [10], [10], [1], [1], [1], [1], [1], [1], [1]],
        [elements] * dims[0]))).T

data = [start_positions]
# define two color
for i in range(data[0].shape[0]):
    rb = np.random.choice([1, 2])
    if rb == 1:
        data[0][i, 4:5] = data[0][i, 4:5] + 1
    if rb == 2:
        data[0][i, 5:6] = data[0][i, 5:6] + 1

particles_list = []
for iteration in range(updates+1):
    # when they have the same coordinate they move togather
    for j in range(data[0].shape[0]):
        xyz = np.random.choice([1, 2, 3])
        minusplus = np.random.choice([-1, 1])
        for k in range(data[0].shape[0]):
            if data[iteration][j, 0:1] == data[iteration][k, 0:1] and data[iteration][j, 1:2] == data[iteration][k, 1:2] and data[iteration][j,2:3] == data[iteration][k, 2:3]:
                data[iteration][j, 6:7] = xyz
                data[iteration][k, 6:7] = data[iteration][j, 6:7]
                data[iteration][j, 7:8] = minusplus
                data[iteration][k, 7:8] = data[iteration][j, 7:8]
            else:
                pass

    previous_data = data[iteration]
    # Movement
    for i in range(data[0].shape[0]):
        xyz = int(data[iteration][i, 6:7])
        minusplus = int(data[iteration][i, 7:8])

        # confine the volume
        box = [0,1,2,3,4,5,6,7,8,9]
        if int(previous_data[i, (xyz-1):xyz] + minusplus) in box:
            previous_data[i, (xyz-1):xyz] = previous_data[i, (xyz-1):xyz] + minusplus
        else:
            previous_data[i, (xyz-1):xyz] = previous_data[i, (xyz-1):xyz]

    # count particles left
    data_array = []
    for i in range(data[0].shape[0]):
        list = (data[iteration][i, 0:3]).tolist() #change array to list
        if list not in data_array:
            data_array.append(list)
    particles = len(data_array)
    data.append(np.array(previous_data))
    particles_list.append(particles)

# let the droplets in the same coordinate become larger
for i in range(updates + 1):
    for j in range(data[0].shape[0]):
        for k in range(data[0].shape[0]):
            if j != k and data[i][j, 0:1] == data[i][k, 0:1] and data[i][j, 1:2] == data[i][k, 1:2] and data[i][j,2:3] == data[i][k, 2:3]:
                data[i][j, 3:4] = data[i][j, 3:4] + 1
            elif j == k:
                pass

# mixing colors
for i in range(updates + 1):
    for j in range(data[0].shape[0]):
        for k in range(data[0].shape[0]):
            if data[i][j, 0:1] == data[i][k, 0:1] and data[i][j, 1:2] == data[i][k, 1:2] and data[i][j,2:3] == data[i][k, 2:3]:
                data[i][j, 8:9] = data[i][j, 8:9] + int(data[i][k, 4:5])
                data[i][j, 9:10] = data[i][j, 9:10] + int(data[i][k, 5:6])
            else:
                pass

#print(data)
# Provide starting angle for the view.
ax.view_init(20, 45)

def scatters(j):
    #clear the previous scatters
    ax.clear()
    # Setting the axes properties
    ax.set_xlim3d([0, 10])
    ax.set_xlabel('X')

    ax.set_ylim3d([0, 10])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0, 10])
    ax.set_zlabel('Z')

    ax.set_title('3D Animated Scatter')
    #update scatter plot

    for i in range(data[0].shape[0]):
        accred = float(data[j][i, 8:9])
        accblue = float(data[j][i, 9:10])
        normred = float(accred / (accred + accblue))
        normblue = float(accblue / (accred + accblue))
        ax.scatter(data[j][i, 0:1], data[j][i, 1:2], data[j][i, 2:3], color=(normred, 0, normblue), s=30 * (data[j][i, 3:4]+1))
    plt.suptitle(r"t = %s" % j + "    particles = %d" % particles_list[j])

def animate(j):
    scatters(j+1)
    return scatters

def init():
    scatters(0)
    return scatters


ani = animation.FuncAnimation(fig=fig, func=animate, frames=updates, init_func=init, interval=6, blit=False,
                              repeat=False)

plt.show()
