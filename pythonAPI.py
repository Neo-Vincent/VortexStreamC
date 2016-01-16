import ctypes
class ns2d:
    ll=ctypes.cdll.LoadLibrary
    lib=ll('vortex.dll') #dll 的文件名
    def __init__(self):
        self.pointer=self.lib.ns2d_new();
    def fdm(self):
        return self.lib.ns2d_fdm(self.pointer)
    def fdm_im(self):
        return self.lib.ns2d_fdm_im(self.pointer)
    def fdm2(self):
        return self.lib.ns2d_fdm2(self.pointer)
    def fvm1_ex(self):
        return self.lib.ns2d_fvm1_ex(self.pointer)
    def fvm2_ex(self):
        return self.lib.ns2d_fvm2_ex(self.pointer)
    def fvm1_im(self):
        return self.lib.ns2d_fvm1_im(self.pointer)
    def fvm2_im(self):
        return self.lib.ns2d_fvm2_im(self.pointer)
    def solve(self,Method=None):
        if Method:
            return self.lib.ns2d_solve1(self.pointer ,ctypes.c_int(Method))
        return self.lib.ns2d_solve0(self.pointer)
    def output(self,filename=None,Dt=None):
        if filename:
            return self.lib.ns2d_output1(self.pointer,ctypes.c_char_p(filename))
        if Dt:
            return self.lib.ns2d_output2(self.pointer,ctypes.c_double(Dt))
        return self.lib.ns2d_output(self.pointer)
    def setDt(self,dt):
        return self.lib.ns2d_setDt(self.pointer, ctypes.c_double(dt))
    def setRe(self,Re):
        return self.lib.ns2d_setRe(self.pointer,ctypes.c_double(Re))
    def setSorwei(self,sorwei):
        return self.lib.ns2d_setSorwei(self.pointer,ctypes.c_double((sorwei)))
    def setDim(Self,a,b):
        return self.lib.ns2d_setDim(self.pointer,ctypes.c_int(a),ctypes,c_int(b))
    def help(self):
        return self.lib.ns2d_help(self.pointer)
    def init(self):
        return self.lib.ns2d_init(self.pointer)
    def setT(self,t):
        return self.lib.ns2d_setT(self.pointer,ctypes.c_double(t))
    def setT0(self,t):
        return self.lib.ns2d_setT0(self.pointer,ctypes.c_double(t))
    def setOutT(self,t):
        return self.lib.ns2d_setOutT(self.pointer,ctypes.c_double(t))
    def isCon(self,isc):
        return self.lib.ns2d_isCon(self.pointer,ctypes.c_bool(isc))
    def setMethod(self,method):
        return self.lib.ns2d_setMethod(self.pointer,ctypes.c_int(method))
    def loadData(self,filename):
        return self.lib.ns2d_loadData(self.pointer,ctypes.c_char_p(filename))
    def setOutName(self,filename):
        return self.lib.ns2d_setOutName(self.pointer,ctypes.c_char_p(filename))
