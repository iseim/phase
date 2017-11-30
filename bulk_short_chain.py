import hoomd
import hoomd.md
import numpy as np


hoomd.context.initialize("");


#initialize molecules' snapshot
num_RNA = 20;
RNA_type = 'CLN3';
if RNA_type == 'CLN3':
	RNA_len = 33;
	protein_len = 8;
	num_protein = np.floor(num_RNA * 21.7);
	num_protein = num_protein.astype(int);
	num_total_particles = RNA_len*num_RNA + protein_len*num_protein;
	cube_side = 22.6 * (num_RNA ** (1./3.));

extent = cube_side;

#particle types: RNA-PRO -> sites on RNA that bind protein
#				 RNA-nb -> neutral sites on RNA
#				 PRO-RNA -> sites on protein that bind RNA
#				 PRO-PRO -> sites on protein that bind protein
#				 PRO-nb -> neutral sites on protein
snapshot = hoomd.data.make_snapshot(N=num_total_particles, box=hoomd.data.boxdim(L=extent),
	particle_types=['RNA-nb', 'RNA-PRO', 'PRO-nb', 'PRO-PRO', 'PRO-RNA'], bond_types=['polymer']);


#create particle positions
from polymer_init_bulk import polymer_init_bulk
pos = polymer_init_bulk(num_RNA, RNA_len, num_protein, protein_len, extent);

snapshot.particles.position[:] = pos;


#define particle types
if RNA_type == 'CLN3':

	#The following code should implemented as a function such that the only input is a vector of the
	#indices of the binding sites on RNA and protein; and the length of the RNA and protein, which
	#are defined earlier in the script as RNA_len and protein_len

	#RNA
	rna_non = np.arange(0, 33);
	rna_sites = np.array([1, 3, 5, 7, 30]);
	rna_non = np.delete(rna_non, rna_sites);

	#protein
	pro_non = np.arange(0, 8);
	pro_sites_pro = np.array([2]);
	pro_sites_RNA = np.array([4, 6]);
	pro_non = np.delete(pro_non, np.array([2, 4, 6]));

	rna_sites2 = rna_sites;
	rna_non2 = rna_non;
	for i in range(num_RNA-1):
		rna_sites2 = np.append(rna_sites2, (i+1)*RNA_len + rna_sites);
		rna_non2 = np.append(rna_non2, (i+1)*RNA_len + rna_non);

	start = num_RNA * RNA_len;
	pro_sites_pro2 = pro_sites_pro;
	pro_sites_RNA2 = pro_sites_RNA;
	pro_non2 = pro_non;
	for i in range(num_protein-1):
		pro_sites_pro2 = np.append(pro_sites_pro2, (i+1)*protein_len + pro_sites_pro);
		pro_sites_RNA2 = np.append(pro_sites_RNA2, (i+1)*protein_len + pro_sites_RNA);
		pro_non2 = np.append(pro_non2, (i+1)*protein_len + pro_non);

	pro_sites_pro2 += start;
	pro_sites_RNA2 += start;
	pro_non2 += start;

snapshot.particles.typeid[rna_non2] = 0;
snapshot.particles.typeid[rna_sites2] = 1;
snapshot.particles.typeid[pro_non2] = 2;
snapshot.particles.typeid[pro_sites_pro2] = 3;
snapshot.particles.typeid[pro_sites_RNA2] = 4;


#no binding sites, just define RNA and protein
#rna = np.arange(0, num_RNA * RNA_len);
#start = num_RNA * RNA_len;
#pro = np.arange(start, start + num_protein * protein_len);
#snapshot.particles.typeid[rna] = 0;
#snapshot.particles.typeid[pro] = 1;


#define bonds
bondnum = num_RNA*(RNA_len-1) + num_protein*(protein_len-1);

snapshot.bonds.resize(bondnum);
bonds = np.zeros(bondnum, dtype='2int64');
index = 0;
for i in range(num_RNA):
	for j in range(RNA_len-1):
		bonds[index][0] = i*RNA_len+j;
		bonds[index][1] = i*RNA_len+j+1;
		index += 1;
start = num_RNA * RNA_len;
for i in range(num_protein):
	for j in range(protein_len-1):
		bonds[index][0] = start + i*protein_len+j;
		bonds[index][1] = start + i*protein_len+j+1;
		index += 1;

snapshot.bonds.group[:] = bonds;


#read snapshot
hoomd.init.read_snapshot(snapshot);


#only interactive force is hard spheres and FENE bond for chains
nl = hoomd.md.nlist.cell();
lj = hoomd.md.pair.lj(r_cut=2**(1./6.), nlist=nl);
lj.pair_coeff.set('RNA-nb', 'RNA-nb', epsilon=1, sigma=1);
lj.pair_coeff.set('RNA-nb', 'RNA-PRO', epsilon=1, sigma=1);
lj.pair_coeff.set('RNA-nb', 'PRO-nb', epsilon=1, sigma=1);
lj.pair_coeff.set('RNA-nb', 'PRO-PRO', epsilon=1, sigma=1);
lj.pair_coeff.set('RNA-nb', 'PRO-RNA', epsilon=1, sigma=1);
lj.pair_coeff.set('RNA-PRO', 'RNA-PRO', epsilon=1, sigma=1);
lj.pair_coeff.set('RNA-PRO', 'PRO-nb', epsilon=1, sigma=1);
lj.pair_coeff.set('RNA-PRO', 'PRO-PRO', epsilon=1, sigma=1);
lj.pair_coeff.set('PRO-nb', 'RNA-PRO', epsilon=1, sigma=1);
lj.pair_coeff.set('PRO-nb', 'PRO-nb', epsilon=1, sigma=1);
lj.pair_coeff.set('PRO-nb', 'PRO-PRO', epsilon=1, sigma=1);
lj.pair_coeff.set('PRO-nb', 'PRO-RNA', epsilon=1, sigma=1);
lj.pair_coeff.set('PRO-PRO', 'PRO-RNA', epsilon=1, sigma=1);
lj.pair_coeff.set('PRO-RNA', 'PRO-PRO', epsilon=1, sigma=1);
lj.pair_coeff.set('PRO-RNA', 'PRO-RNA', epsilon=1, sigma=1);

lj.pair_coeff.set('RNA-PRO', 'PRO-RNA', epsilon=0, sigma=1);
lj.pair_coeff.set('PRO-PRO', 'PRO-PRO', epsilon=0, sigma=1);


#add in attractive forces
morse = hoomd.md.pair.morse(r_cut=5, nlist=nl);
morse.pair_coeff.set('RNA-PRO', 'PRO-RNA', D0=14, alpha=3, r0=1);
morse.pair_coeff.set('PRO-PRO', 'PRO-PRO', D0=7, alpha=3, r0=1);

morse.pair_coeff.set('RNA-nb', 'RNA-nb', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('RNA-nb', 'RNA-PRO', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('RNA-nb', 'PRO-nb', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('RNA-nb', 'PRO-PRO', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('RNA-nb', 'PRO-RNA', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('RNA-PRO', 'RNA-PRO', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('RNA-PRO', 'PRO-nb', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('RNA-PRO', 'PRO-PRO', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('PRO-nb', 'RNA-PRO', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('PRO-nb', 'PRO-nb', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('PRO-nb', 'PRO-PRO', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('PRO-nb', 'PRO-RNA', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('PRO-PRO', 'PRO-RNA', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('PRO-RNA', 'PRO-PRO', D0=0, alpha=3, r0=1);
morse.pair_coeff.set('PRO-RNA', 'PRO-RNA', D0=0, alpha=3, r0=1);


#define FENE bond
#fene = hoomd.md.bond.fene(name="fenebond");
#fene.bond_coeff.set('polymer', k=30, r0=1.5, sigma=1, epsilon=1);
#fene.disable()

#define harmonic bond
harm = hoomd.md.bond.harmonic(name="harmbond");
harm.bond_coeff.set('polymer', k=400, r0=2**(1./6.));


#select integrator
integrator_mode = hoomd.md.integrate.mode_standard(dt=.00025);
br = hoomd.md.integrate.brownian(group=hoomd.group.all(), kT=1, seed=1);


#write output
hoomd.dump.gsd('dilute.gsd', period=500, group=hoomd.group.all(), overwrite=True);
#logger = hoomd.analyze.log(filename="thermo2.log", period=1000, quantities=["temperature"], overwrite=True)



hoomd.run(10e6)









