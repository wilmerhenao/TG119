
import numpy as np
import scipy.io as io
import scipy.sparse as sps
from gurobipy import *


class fmo_data(object):
    def __init__(self,inputFilename):
        #Reads in map file
        matFile = io.loadmat(inputFilename)
        self.nVox, self.nBeams, self.nBPB = int(matFile['nVox']), int(matFile['nBeams']), np.array(matFile['nBPB'].flatten())
        self.nBix, structures, self.nStruct = self.nBPB.sum(), matFile['structures'], int(matFile['nStruct'])
        self.gantry, self.couch, self.dataFolder = matFile['ga'].flatten(), matFile['ca'].flatten(), str(matFile['dataFolder'])
        self.structs = [fmo_struct(structures[0,b][0]) for b in range(self.nStruct)]
        dHolder = []
        for b in range(self.nBeams):
            beamDFile = io.loadmat(self.dataFolder + 'Gantry' + str(self.gantry[b])+'_Couch' + str(self.couch[b]) + '_Dnew.mat')
            dHolder.append(beamDFile['Dnew'])
        self.D = sps.hstack([dHolder[d] for d in range(len(dHolder))]).tocsc()
        self.printAllData()
        
    def printAllData(s):
        #Prints data
        print('voxels=',s.nVox,'beams=',s.nBeams,'beamlets=',s.nBPB.sum())
        print('gantry=',s.gantry,'couch=', s.couch)
        
class fmo_struct(object):
    def __init__(self,dataStruct):
        self.name = str(dataStruct[0][0])
        self.target = float(dataStruct[1][0])
        self.meanUB = float(dataStruct[2][0])
        self.doseLB = float(dataStruct[3][0])
        self.doseUB = float(dataStruct[4][0])
        self.voxels = dataStruct[5]-1




class fmo_model(object):
    def __init__(self,data):

        print('initializing model')
        self.m = Model('FMO')

        print('generating z variables')
        self.z = [self.m.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in xrange(data.nVox)]
        self.x = []
        self.m.update()


        print('initializing dose constraint')
        self.doseConstr = [self.m.addConstr(-self.z[j],GRB.EQUAL,0) for j in xrange(data.nVox)]
        self.m.update()

        print( 'generating x variables and populating dose constraint')
        for i in range(data.nBix):
            xCol = Column()
            indices,dcol = data.D[:,i].indices, data.D[:,i].data
            for d in range(len(indices)):
                xCol.addTerms(dcol[d],self.doseConstr[indices[d]-1])
            self.x.append(self.m.addVar(vtype=GRB.CONTINUOUS,column = xCol))
        self.m.update()


        for s in range(len(data.structs)):
            print('building constraints for structure',data.structs[s].name)
            ##build mean bound
            if data.structs[s].meanUB>0:
                data.structs[s].meanVar = self.m.addVar(vtype=GRB.CONTINUOUS); self.m.update()
                self.m.addConstr(len(data.structs[s].voxels)*data.structs[s].meanVar,GRB.EQUAL,quicksum([self.z[j] for j in data.structs[s].voxels]))
                self.m.update()

            ##build upper dose bound
            if data.structs[s].doseUB>0:
                for j in data.structs[s].voxels:
                    self.z[j].setAttr('UB',data.structs[s].doseUB)
                self.m.update()
            ##build lower dose bound
            if data.structs[s].doseLB>0:
                for j in data.structs[s].voxels:
                    self.z[j].setAttr('LB',data.structs[s].doseLB)
                self.m.update()

        print('building objective (maximize mean dose to ptv)')
        self.objExpr = LinExpr()
        for s in range(len(data.structs)):
            if data.structs[s].target==1:
                if not hasattr(data.structs[s],'meanVar'):
                    data.structs[s].meanVar = self.m.addVar(vtype=GRB.CONTINUOUS); self.m.update()
                    self.m.addConstr(len(data.structs[s].voxels)*data.structs[s].meanVar,GRB.EQUAL,quicksum([self.z[j] for j in data.structs[s].voxels]))
                    self.m.update()
                self.objExpr += data.structs[s].meanVar
        self.m.setObjective(self.objExpr, GRB.MAXIMIZE)
        self.m.update()
        
    def solveandoutput(self,data):
        self.m.optimize()
        self.writeSolution(data,'intensitiesOut.mat')

        #print 'writing out model';self.m.write('model.lp')
    def writeSolution(self,data,outName):
        #This saves the apertures to a matlab cell array
        io.savemat(outName,{'x':[self.x[i].X for i in range(len(self.x)) ],'subdose':[self.z[j].X for j in range(len(self.z))]})

data = fmo_data('tg119pythondata.mat')

model = fmo_model(data)

model.solveandoutput(data)
