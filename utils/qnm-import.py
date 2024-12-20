import numpy as np 
import qnm
import pickle

#Get valid m values for a given l_max
def m_from_l(lmax):
    return [i for i in range(-lmax,lmax+1)]

#Load the data
qnm.download_data()
data_dir = qnm.cached.get_cachedir() / 'data'

s_vals =["-2","-1"]
l_vals = [str(i) for i in range(8)]
n_vals = [str(i) for i in range(6)]

output = 0

for l in l_vals:
    print(f"This is l={l}")
    m_vals = m_from_l(int(l))
    for m in m_vals:
        for n in n_vals:
            for s in s_vals:
                file_loc = f"{f_loc}s{s}_l{l}_m{m}_n{n}.pickle"
                with open(file_loc ,"rb") as f:
                    x = pickle.load(f)
                    out = [["a","re_omega","im_omega"]]
                    for a,omega in zip(x.a,x.omega):
                        out.append([a,omega.real,omega.imag])
                    if output:
                        np.savetxt(f"data_as_txt/s{s}_l{l}_m{m}_n{n}.txt",np.array(out),fmt="%s",delimiter=",")
