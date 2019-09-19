import matplotlib.pyplot as plt
import csv
import sys

ticks = []
tSeries = []
mProfile_ff = []
mProfile = []
error = []

with open('../configs/' +  sys.argv[1] + '.cfg','r') as f:
    #w, h = [int(x) for x in next(f).split()] # read first line
    array = []
    for line in range(4): # read rest of lines
        array.append([int(x) for x in next(f).split()])

with open('../results/result_' +  sys.argv[1] + '.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        ticks.append(int(row[0]))
        tSeries.append(float(row[1]))
        mProfile_ff.append(float(row[2]))
        mProfile.append(float(row[4]))
        if float(row[6]) > 100:
            error.append(float(0))
        else:
            error.append(float(row[6]))


plt.figure(figsize=(11,8))

plt.subplot(2, 2, 1)
plt.plot(ticks, tSeries, color='b')
plt.title('Original time series')

plt.subplot(2, 2, 4)
plt.plot(ticks, mProfile_ff, color='g')
plt.title('SCRIMP FlexFloat')

plt.subplot(2, 2, 2)
plt.plot(ticks, mProfile, color='g')
plt.title('SCRIMP double precision')

plt.subplot(2, 2, 3)
plt.plot(ticks, error, color='r')
plt.title('Error %')
plt.ylim((0, 100))


plt.text(0.13, 0.03, 'FF parameters [exp, man] => distance=' + str(array[0]) + '; dotprod=' +str(array[1]) + '; stats='+ str(array[2]) + '; profile='+ str(array[3]), fontsize=12, transform=plt.gcf().transFigure)

plt.legend()
plt.show()
