import numpy as np
import matplotlib.pyplot as plt

def read_grid(file_input,is_curv=False):
    # Reading of the grid
    # -------------------
    print("Reading grid...")

    f = open(file_input,'r')
    arg = np.fromfile(f, dtype=('>i4'), count=1)
    nx  = np.fromfile(f, dtype=('>i4'), count=1)[0]
    #===============================================================================
    arg = np.fromfile(f, dtype=('>i8'), count=1)
    ny  = np.fromfile(f, dtype=('>i4'), count=1)[0]
    #===============================================================================
    arg = np.fromfile(f, dtype=('>i8'), count=1)
    nz  = np.fromfile(f, dtype=('>i4'), count=1)[0]
    #===============================================================================
    arg = np.fromfile(f, dtype=('>i8'), count=1)
    if is_curv:
        x = np.zeros((nx,ny),dtype='float',order='F')
        for j in range(ny):
            x[:,j] = np.fromfile(f, dtype=('>f8'), count=nx)
        #===============================================================================
        arg = np.fromfile(f, dtype=('>i8'), count=1)
        y = np.zeros((nx,ny),dtype='float',order='F')
        for j in range(ny):
            y[:,j] = np.fromfile(f, dtype=('>f8'), count=nx)
        #===============================================================================
        arg = np.fromfile(f, dtype=('>i8'), count=1)
        z = np.fromfile(f, dtype=('>f8'), count=ny)
        #===============================================================================
    else:
        x = np.zeros((nx),dtype='float',order='F')
        x = np.fromfile(f, dtype=('>f8'), count=nx)
        #===============================================================================
        arg = np.fromfile(f, dtype=('>i8'), count=1)
        y = np.zeros((ny),dtype='float',order='F')
        y = np.fromfile(f, dtype=('>f8'), count=ny)
        #===============================================================================
        arg = np.fromfile(f, dtype=('>i8'), count=1)
        z = np.zeros((ny),dtype='float',order='F')
        z = np.fromfile(f, dtype=('>f8'), count=nz)
        #===============================================================================

    return nx,ny,nz,x,y,z

def read_restart(file_input,nx,ny,nz,is_2D=False):
    print("Reading restart...")
    f=open(file_input,'r')

    ro  = np.fromfile(f, dtype=('<f8'), count=nx*ny*nz).reshape((nx,ny,nz), order='F')
    rou = np.fromfile(f, dtype=('<f8'), count=nx*ny*nz).reshape((nx,ny,nz), order='F')
    rov = np.fromfile(f, dtype=('<f8'), count=nx*ny*nz).reshape((nx,ny,nz), order='F')
    row = np.fromfile(f, dtype=('<f8'), count=nx*ny*nz).reshape((nx,ny,nz), order='F')
    roe = np.fromfile(f, dtype=('<f8'), count=nx*ny*nz).reshape((nx,ny,nz), order='F')

    if is_2D:
        ro_,rou_,rov_,row_,roe_ = ro,rou,rov,row,roe
        ro,rou,rov,row,roe = np.ndarray((nx,ny)),np.ndarray((nx,ny)),np.ndarray((nx,ny)),np.ndarray((nx,ny)),np.ndarray((nx,ny))
        for j in range(ny):
            for i in range(nx):
                ro[i][j]  = ro_[i][j][0]
                rou[i][j] = rou_[i][j][0]
                rov[i][j] = rov_[i][j][0]
                row[i][j] = row_[i][j][0]
                roe[i][j] = roe_[i][j][0]
        ro = np.transpose(ro)
        rou = np.transpose(rou)
        rov = np.transpose(rov)
        row = np.transpose(row)
        roe = np.transpose(roe)

    return ro,rou,rov,row,roe


def read_stats_hit(file_input):
    print("Reading stats...")
    stats = {'Tstar':[],'rho':[],'u':[],'v':[],'w':[],'p':[],'T':[],'e':[],'h':[],'c':[],\
             's':[],'Mt':[],'0.5*q':[],'g':[],'mu':[],'cok':[],'cp':[],'cv':[],'pr':[],\
             'eck':[],'div':[],'rho*dux':[],'rho*duy':[],'rho*duz':[],'rho*dvx':[],\
             'rho*dvy':[],'rho*dvz':[],'rho*dwx':[],'rho*dwy':[],'rho*dwz':[],'p*div':[],\
             'rho*div':[],'b1':[],'b2':[],'b3':[],'Ttot':[],'vide':[],'rhou':[],'rhov':[],\
             'rhow':[],'rhoe':[],'rho*T':[],'rho^2':[],'u^2':[],'v^2':[],'w^2':[],'uv':[],\
             'uw':[],'vw':[],'vT':[],'p^2':[],'T^2':[],'e^2':[],'h^2':[],'c^2':[],'s^2':[],\
             'Mt^2':[],'g^2':[],'mu^2':[],'cok^2':[],'cv^2':[],'cp^2':[],'pr^2':[],'eck^2':[],\
             'p*u':[],'p*v':[],'T*u':[],'T*v':[],'s*u':[],'s*v':[],'p*rho':[],'T*rho':[],'h*rho':[],\
             'T*p':[],'p*s':[],'T*s':[],'rho*s':[],'g*rho':[],'g*p':[],'g*s':[],'g*T':[],'g*u':[],\
             'g*v':[],'p*dux':[],'p*dvy':[],'p*dwz':[],'p*duy':[],'p*dvx':[],'div^2':[],'rho*div^2':[],\
             'dux^2':[],'duy^2':[],'duz^2':[],'dvx^2':[],'dvy^2':[],'dvz^2':[],'dwx^2':[],'dwy^2':[],\
             'dwz^2':[],'b1^2':[],'b2^2':[],'b3^2':[],'rho*b1':[],'rho*b2':[],'rho*b3':[],'Ttot^2':[],\
             'vide2':[],'rho*u^2':[],'rho*v^2':[],'rho*w^2':[],'rho*T^2':[],'rho*b1^2':[],'rho*b2^2':[],\
             'rho*b3^2':[],'rho*vw':[],'rho*vT':[],'rho*u^2*v':[],'rho*v^3':[],'rho*w^2*v':[],'rho*v^2*u':[],\
             'rho*dux^2':[],'rho*dvy^2':[],'rho*dwz^2':[],'rho*duy*dvx':[],'rho*duz*dwx':[],'rho*dvz*dwy':[],\
             'u^3':[],'p^3':[],'u^4':[],'p^4':[],'Frhou':[],'Frhov':[],'Frhow':[],'Grhov':[],'Grhow':[],\
             'Hrhow':[],'Frhovu':[],'Frhovv':[],'Grhovu':[],'Grhovv':[],'Grhoww':[],'Frhou_dux':[],\
             'Frhou_dvx':[],'Frhov_dux':[],'Frhov_duy':[],'Frhov_dvx':[],'Frhov_dvy':[],'Frhow_dwx':[],\
             'Grhov_duy':[],'Grhov_dvy':[],'Grhow_duz':[],'Grhow_dvz':[],'Grhow_dwy':[],'Hrhow_dwz':[],\
             'cok*dTx':[],'cok*dTy':[],'cok*dTz':[],'h*u':[],'h*v':[],'h*w':[],'rho*h*u':[],\
             'rho*h*v':[],'rho*h*w':[],'rho*u^3':[],'rho*v^3':[],'rho*w^3':[],'vide3':[],'dux^3':[],'duy^3':[],\
             'duz^3':[],'dvx^3':[],'dvy^3':[],'dvz^3':[],'dwx^3':[],'dwy^3':[],'dwz^3':[],'dux^4':[],'duy^4':[],\
             'duz^4':[],'dvx^4':[],'dvy^4':[],'dvz^4':[],'dwx^4':[],'dwy^4':[],'dwz^4':[],'mu_rho':[],'4_3*mu*div^2':[],\
             'mu*vrt^2':[],'maxg_ming':[],'rho*uv':[],'Frhow_duz':[],'Frhow_dvz':[]}

    # Reading of the stats
    # --------------------
    f = open(file_input,'r')
    for line in f.readlines():
            save = line.split()
            stats['Tstar'].append(float(save[0]))
            stats['rho'].append(float(save[1]))
            stats['u'].append(float(save[2]))
            stats['v'].append(float(save[3]))
            stats['w'].append(float(save[4]))
            stats['p'].append(float(save[5]))
            stats['T'].append(float(save[6]))
            stats['e'].append(float(save[7]))
            stats['h'].append(float(save[8]))
            stats['c'].append(float(save[9]))
            stats['s'].append(float(save[10]))
            stats['Mt'].append(float(save[11]))
            stats['0.5*q'].append(float(save[12]))
            stats['g'].append(float(save[13]))
            stats['mu'].append(float(save[14]))
            stats['cok'].append(float(save[15]))
            stats['cp'].append(float(save[16]))
            stats['cv'].append(float(save[17]))
            stats['pr'].append(float(save[18]))
            stats['eck'].append(float(save[19]))
            stats['div'].append(float(save[20]))
            stats['rho*dux'].append(float(save[21]))
            stats['rho*duy'].append(float(save[22]))
            stats['rho*duz'].append(float(save[23]))
            stats['rho*dvx'].append(float(save[24]))
            stats['rho*dvy'].append(float(save[25]))
            stats['rho*dvz'].append(float(save[26]))
            stats['rho*dwx'].append(float(save[27]))
            stats['rho*dwy'].append(float(save[28]))
            stats['rho*dwz'].append(float(save[29]))
            stats['p*div'].append(float(save[30]))
            stats['rho*div'].append(float(save[31]))
            stats['b1'].append(float(save[32]))
            stats['b2'].append(float(save[33]))
            stats['b3'].append(float(save[34]))
            stats['Ttot'].append(float(save[35]))
            stats['vide'].append(float(save[36]))
            stats['rhou'].append(float(save[37]))
            stats['rhov'].append(float(save[38]))
            stats['rhow'].append(float(save[39]))
            stats['rhoe'].append(float(save[40]))
            stats['rho*T'].append(float(save[41]))
            stats['rho^2'].append(float(save[42]))
            stats['u^2'].append(float(save[43]))
            stats['v^2'].append(float(save[44]))
            stats['w^2'].append(float(save[45]))
            stats['uv'].append(float(save[46]))
            stats['uw'].append(float(save[47]))
            stats['vw'].append(float(save[48]))
            stats['vT'].append(float(save[49]))
            stats['p^2'].append(float(save[50]))
            stats['T^2'].append(float(save[51]))
            stats['e^2'].append(float(save[52]))
            stats['h^2'].append(float(save[53]))
            stats['c^2'].append(float(save[54]))
            stats['s^2'].append(float(save[55]))
            stats['Mt^2'].append(float(save[56]))
            stats['g^2'].append(float(save[57]))
            stats['mu^2'].append(float(save[58]))
            stats['cok^2'].append(float(save[59]))
            stats['cv^2'].append(float(save[60]))
            stats['cp^2'].append(float(save[61]))
            stats['pr^2'].append(float(save[62]))
            stats['eck^2'].append(float(save[63]))
            stats['p*u'].append(float(save[64]))
            stats['p*v'].append(float(save[65]))
            stats['T*u'].append(float(save[66]))
            stats['T*v'].append(float(save[67]))
            stats['s*u'].append(float(save[68]))
            stats['s*v'].append(float(save[69]))
            stats['p*rho'].append(float(save[70]))
            stats['T*rho'].append(float(save[71]))
            stats['h*rho'].append(float(save[72]))
            stats['T*p'].append(float(save[73]))
            stats['p*s'].append(float(save[74]))
            stats['T*s'].append(float(save[75]))
            stats['rho*s'].append(float(save[76]))
            stats['g*rho'].append(float(save[77]))
            stats['g*p'].append(float(save[78]))
            stats['g*s'].append(float(save[79]))
            stats['g*T'].append(float(save[80]))
            stats['g*u'].append(float(save[81]))
            stats['g*v'].append(float(save[82]))
            stats['p*dux'].append(float(save[83]))
            stats['p*dvy'].append(float(save[84]))
            stats['p*dwz'].append(float(save[85]))
            stats['p*duy'].append(float(save[86]))
            stats['p*dvx'].append(float(save[87]))
            stats['div^2'].append(float(save[88]))
            stats['rho*div^2'].append(float(save[89]))
            stats['dux^2'].append(float(save[90]))
            stats['duy^2'].append(float(save[91]))
            stats['duz^2'].append(float(save[92]))
            stats['dvx^2'].append(float(save[93]))
            stats['dvy^2'].append(float(save[94]))
            stats['dvz^2'].append(float(save[95]))
            stats['dwx^2'].append(float(save[96]))
            stats['dwy^2'].append(float(save[97]))
            stats['dwz^2'].append(float(save[98]))
            stats['b1^2'].append(float(save[99]))
            stats['b2^2'].append(float(save[100]))
            stats['b3^2'].append(float(save[101]))
            stats['rho*b1'].append(float(save[102]))
            stats['rho*b2'].append(float(save[103]))
            stats['rho*b3'].append(float(save[104]))
            stats['Ttot^2'].append(float(save[105]))
            stats['vide2'].append(float(save[106]))
            stats['rho*u^2'].append(float(save[107]))
            stats['rho*v^2'].append(float(save[108]))
            stats['rho*w^2'].append(float(save[109]))
            stats['rho*T^2'].append(float(save[110]))
            stats['rho*b1^2'].append(float(save[111]))
            stats['rho*b2^2'].append(float(save[112]))
            stats['rho*b3^2'].append(float(save[113]))
            stats['rho*uv'].append(float(save[114]))
            stats['rho*vw'].append(float(save[115]))
            stats['rho*vT'].append(float(save[116]))
            stats['rho*u^2*v'].append(float(save[117]))
            stats['rho*v^3'].append(float(save[118]))
            stats['rho*w^2*v'].append(float(save[119]))
            stats['rho*v^2*u'].append(float(save[120]))
            stats['rho*dux^2'].append(float(save[121]))
            stats['rho*dvy^2'].append(float(save[122]))
            stats['rho*dwz^2'].append(float(save[123]))
            stats['rho*duy*dvx'].append(float(save[124]))
            stats['rho*duz*dwx'].append(float(save[125]))
            stats['rho*dvz*dwy'].append(float(save[126]))
            stats['u^3'].append(float(save[127]))
            stats['p^3'].append(float(save[128]))
            stats['u^4'].append(float(save[129]))
            stats['p^4'].append(float(save[130]))
            stats['Frhou'].append(float(save[131]))
            stats['Frhov'].append(float(save[132]))
            stats['Frhow'].append(float(save[133]))
            stats['Grhov'].append(float(save[134]))
            stats['Grhow'].append(float(save[135]))
            stats['Hrhow'].append(float(save[136]))
            stats['Frhovu'].append(float(save[137]))
            stats['Frhovv'].append(float(save[138]))
            stats['Grhovu'].append(float(save[139]))
            stats['Grhovv'].append(float(save[140]))
            stats['Grhoww'].append(float(save[141]))
            stats['Frhou_dux'].append(float(save[142]))
            stats['Frhou_dvx'].append(float(save[143]))
            stats['Frhov_dux'].append(float(save[144]))
            stats['Frhov_duy'].append(float(save[145]))
            stats['Frhov_dvx'].append(float(save[146]))
            stats['Frhov_dvy'].append(float(save[147]))
            stats['Frhow_duz'].append(float(save[148]))
            stats['Frhow_dvz'].append(float(save[149]))
            stats['Frhow_dwx'].append(float(save[150]))
            stats['Grhov_duy'].append(float(save[151]))
            stats['Grhov_dvy'].append(float(save[152]))
            stats['Grhow_duz'].append(float(save[153]))
            stats['Grhow_dvz'].append(float(save[154]))
            stats['Grhow_dwy'].append(float(save[155]))
            stats['Hrhow_dwz'].append(float(save[156]))
            stats['cok*dTx'].append(float(save[157]))
            stats['cok*dTy'].append(float(save[158]))
            stats['cok*dTz'].append(float(save[159]))
            stats['h*u'].append(float(save[160]))
            stats['h*v'].append(float(save[161]))
            stats['h*w'].append(float(save[162]))
            stats['rho*h*u'].append(float(save[163]))
            stats['rho*h*v'].append(float(save[164]))
            stats['rho*h*w'].append(float(save[165]))
            stats['rho*u^3'].append(float(save[166]))
            stats['rho*v^3'].append(float(save[167]))
            stats['rho*w^3'].append(float(save[168]))
            stats['vide3'].append(float(save[169]))
            stats['dux^3'].append(float(save[170]))
            stats['duy^3'].append(float(save[171]))
            stats['duz^3'].append(float(save[172]))
            stats['dvx^3'].append(float(save[173]))
            stats['dvy^3'].append(float(save[174]))
            stats['dvz^3'].append(float(save[175]))
            stats['dwx^3'].append(float(save[176]))
            stats['dwy^3'].append(float(save[177]))
            stats['dwz^3'].append(float(save[178]))
            stats['dux^4'].append(float(save[179]))
            stats['duy^4'].append(float(save[180]))
            stats['duz^4'].append(float(save[181]))
            stats['dvx^4'].append(float(save[182]))
            stats['dvy^4'].append(float(save[183]))
            stats['dvz^4'].append(float(save[184]))
            stats['dwx^4'].append(float(save[185]))
            stats['dwy^4'].append(float(save[186]))
            stats['dwz^4'].append(float(save[187]))
            stats['mu_rho'].append(float(save[188]))
            stats['4_3*mu*div^2'].append(float(save[189]))
            stats['mu*vrt^2'].append(float(save[190]))
            stats['maxg_ming'].append(float(save[191]))
    return stats
