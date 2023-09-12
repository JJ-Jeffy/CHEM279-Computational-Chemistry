import matplotlib.pyplot as plt


def calculate_slope(x, y):
    # calculate the slope of the line
    slope = (y[1] - y[0])/(x[1] - x[0])
    return slope

def plot_data(x, y, title, x_label, y_label):
    plt.plot(x, y)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend(['Forward', 'Central'])
    # calculate the slope of the line
    slope = calculate_slope(x, y)
    print(slope)
    # save the file as a png
    plt.savefig('Error.png')


def read_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    # skip the first line 
    lines = lines[1:]
    x = []
    y1 = []
    y2 = []
    for line in lines:
        line = line.strip()
        if line:
            line = line.split(' ')
            x.append(float(line[0]))
            y1.append(float(line[1]))
            y2.append(float(line[2]))
    return x, y1, y2 

x, y1, y2 = read_file('error.txt')
plot_data(x, y1, 'Error_Forward', 'log(h)', 'log(Error)')
plot_data(x, y2, 'Error_central', 'log(h)', 'log(Error)')