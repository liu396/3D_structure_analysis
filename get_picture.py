import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('45RD_matrix.pdf')

# plt.imshow(np.random.random((50,50)));
# plt.colorbar()
# plt.show()

#count=np.array([1,2,3,4]);
#grid=np.reshape(count,(2,2));
#plt.imshow(grid,interpolation='none',cmap=plt.get_cmap('jet'))
#plt.colorbar(orientation='vertical')
#plt.show()
#exit();

num_rows=11;

k1=[]
k2=[]
count=[]
for i in range(num_rows*num_rows):
    count.append(0)
    
print(len(count))
    
with open('x.vs.plane.txt') as f:
    print("get in")
    f.readline()
    #f.readline()
    for line in f:
        i, j, p = [float(x) for x in line.split()] # read first line
        i=int(i)
        j=int(j)
        if (i>j and p!=0): 
            print(i,j,p)
        k1.append(i)
        k2.append(j)
        #print(i*num_rows+num_rows-j-1)
        count[(num_rows-j-1)*num_rows+i]=p
        
k1=np.array(k1)
k2=np.array(k2)
count=np.array(count)
print(len(k1))
print(len(k2))
print(len(count))
grid=count.reshape(num_rows,num_rows)
plt.imshow(grid, extent=(0,1.1,0,1.1),interpolation='none',cmap=plt.get_cmap('jet'))
#plt.imshow(grid,interpolation='bilinear',cmap=plt.get_cmap('inferno'))
#plt.matshow(grid,cmap=plt.get_cmap('inferno'))
plt.xlabel('X',fontsize=15)
plt.ylabel('Z',fontsize=15)
plt.colorbar(orientation='vertical')
#plt.clim(0,0.6)
plt.savefig(pp, format='pdf')
pp.close()
plt.show()
