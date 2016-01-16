import matplotlib.pyplot as plt
import numpy as np
import time
class fluentData():
    filename='fluent400400.dat'

    def fileDataRead(self,filename=None):
        data=[]
        if filename:
            self.filename=filename
        f=open(self.filename)
        newfile=self.filename+'_plot'
        g=open(newfile,'w')

        line=f.readlines()
        self.line=line
        bloc=line[0].replace('\n','').split(',')
        a='variables="'
        for i in range(len(bloc)-2):
            a+=bloc[i+1]+'","'
        a=a+bloc[-1]+'"\n'
        #N=str(np.sqrt(float(line[-1])))
        N='101'
        b='zone i='+N+'j='+N+' f=point\n'
        g.write(a)
        g.write(b)

    def readData(self,filename=None):
        if filename:
            self.filename=filename
        f=open(self.filename)
        lines=f.readlines()
        f.close()
        newfile=self.filename+'.plt'
        bloc=lines[0].replace('\n','').split(',')
        a='variables="'
        for i in range(len(bloc)-2):
            a+=bloc[i+1]+'","'
        a=a+bloc[-1]+'"\n'
        a=a.replace('x-coordinate','x')\
            .replace('y-coordinate','y')\
            .replace('pressure-coefficient','p')\
            .replace("x-velocity",'u')\
            .replace("y-velocity",'v')\
            .replace('stream-function','stream line')
        data=[]
        setx=set()
        sety=set()
        print("loading data\n")
        for i in range(1,len(lines)):
            block=lines[i].replace('\n','').split(',')
            if len(block)!=0:
                values=list(map(lambda x: float(x),block[1:]))
                if values not in data:
                    data.append(values)
                setx.add(values[0])
                sety.add(values[1])
        nx=len(setx)
        ny=len(sety)
        N=len(data)
        #check the count
        if N!=nx*ny:
            print("counting nx ny error!\n\
nx={0},ny={1},N={2}\n\
but the file will be generated anyway!\n".format(nx,ny,N))
        # rearrange the line from x=min(x) to max(x)
        print("rearrange data\n")
        t0=time.time()
        print("t0=",time.time());
        data.sort()
        print("t1-t0=",time.time()-t0)
        a+='zone i='+str(nx)+'  j='+str(ny)+'  f=point\n'
        for i in data:
            a+=','.join(list(map(lambda x:str(x),i)))+'\n'
        g=open(newfile,mode='w')
        g.write(a)
        g.close()
'''
        for line in content:

            if line[0]!='#':



                line.replace('\n','')
                blok=line.split(',')
                dim=len(blok)
                vector=[]
                for xy in blok:
                    try:
                        vector.append(float(xy))
                    except ValueError:
                        pass
                if len(vector)!=0:
                    data.append(vector)
        return data

def draw(data):
    # x=0.5 u as a function of y
    Fig1=plt.figure(1)
    nb=len(data)
    x=[]
    y=[]
    plt.title('x = 0.5 u as a function of y')# give plot a title
    plt.xlabel('u at x=0.5')# make axis labels
    plt.ylabel('y')

    for i in range(nb):
        x.append(data[i][2])
        y.append(data[i][1])
    plt.plot(x,y,'ro-')
    Fig1.show()
    Fig1.savefig("x = 0.5 u as a function of y.png")
    
    Fig2=plt.figure(2)
    nb=len(data)
    x=[]
    y=[]
    plt.title('y = 0.5 v as a function of x')# give plot a title
    plt.xlabel('x')# make axis labels
    plt.ylabel('v at y=0.5')

    for i in range(nb):
        x.append(data[i][1])
        y.append(data[i][3])
    plt.plot(x,y,'ro-')
    Fig2.show()
    Fig2.savefig("y = 0.5 v as a function of x'.png")



a='X=0.5_uvp.dat'
data=fileDataRead(a)
draw(data)

 variables= "  x  "  , "  y","u","v","vorticity","stream line","P"
  zone i=101   j=101 f=point 

'''
a=fluentData()
a.readData()
