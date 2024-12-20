import numpy as np 
import qnm
import pickle
import h5py

#Get valid m values for a given l_max
def m_from_l(lmax):
    return [i for i in range(-lmax,lmax+1)]

#Load the data
qnm.download_data()
data_dir = qnm.cached.get_cachedir() / "data"

output = 0

for s in [-2, -1]:
  with h5py.File(f"QNM_s{s}.h5", 'w') as qnm_h5:
    qnm_h5.attrs["s"] = s

    for n in range(0,7):
      for l in range(abs(s),8):
        for m in range(-l,l+1):
          print("Processing: ", s, n, l, m)
          try:
            with open(data_dir / f"s{s}_l{l}_m{m}_n{n}.pickle", "rb") as f:
              qnmdata = pickle.load(f)
              qnm_h5.create_dataset(f"l{l}/m{m}/n{n}/omega", data=qnmdata.omega, shuffle=True, compression='gzip', compression_opts=6, track_times=False)
              qnm_h5.create_dataset(f"l{l}/m{m}/n{n}/a", data=qnmdata.a, shuffle=True, compression='gzip', compression_opts=6, track_times=False)
          except:
            print("Failed to create dataset, skipping.")
