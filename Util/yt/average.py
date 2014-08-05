import yt
import pylab

ds = yt.load("plt164582")

dd = ds.all_data()

p = yt.create_profile(dd, "z", "temperature", 
                      n_bins=ds.domain_dimensions[2],
                      weight_field="cell_mass",
                      logs={"z": False})


#pylab.plot(p.x, p[("gas","temperature")])
pylab.plot(p.x, p.variance[("gas","temperature")]/p[("gas","temperature")])

ax = pylab.gca()
ax.set_yscale("log")

pylab.savefig("average.png")


