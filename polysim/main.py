import hoomd
import hoomd.md
import numpy as np
import gsd.hoomd
import random
import fresnel
import PIL
import IPython.display

def simulate_polymer_dpd(N, A, gamma = 1.0, sigma = 1.0, dt=0.01, steps = 1e5, kT = 0.8, trajectory_filename = 'trajectory.gsd', traj_frames = 5, ):

    # N is the number of monomers per chain (i.e. degree of polymerization)
    # lp is the persistence length of the polymer 
    # sigma is the diameter of bead 
    # dt is the step size
    # steps is the number of timesteps
    # kT is the simulation temperature

    ## Initialize sytem
    hoomd.util.quiet_status()
    hoomd.context.initialize("");
    hoomd.option.set_notice_level(0);  

    ## Set dimensions of simulation box to accomodate size of polymer
    L = N*sigma # box length (check units)

    # Build simulation snapshot
    snapshot = hoomd.data.make_snapshot(N = N,
                                        box=hoomd.data.boxdim(Lx=L, Ly=L, Lz=L),
                                        particle_types=['A'],
                                        bond_types=['polymer']);

    ## Build linear chain coordiantes
    positions = []
    for i in np.arange(N):
        x = -L/2 + 1/2 + i
        y = 0
        z = 0
        positions.append([x,y,z])

    snapshot.particles.position[:] = positions

    ## Define particle type 
    snapshot.particles.typeid[0:N]=0;

    ## Define particle bonds 
    snapshot.bonds.resize(N - 1);
    bond_pairs = [[i, i+1] for i in np.arange(N-1)]
    snapshot.bonds.group[:] = bond_pairs

    # Check build
    hoomd.init.read_snapshot(snapshot);

    # Define Force-Field
    ## Find Neighboring Particles
    nl = hoomd.md.nlist.cell();

    ## Set DPD cut-off
    dpd = hoomd.md.pair.dpd(r_cut=3.5*sigma, nlist=nl, kT=kT, seed=random.randint(1,10**4));

    ## Set DPD potential
    dpd.pair_coeff.set('A', 'A', A=A, gamma = gamma);

    ## Set DPD exclusions
    nl.reset_exclusions(exclusions = []);

    ## Set Harmonic Bond Potential
    harmonic = hoomd.md.bond.harmonic();
    harmonic.bond_coeff.set('polymer', k=100.0, r0=sigma);

    ## Select integrator
    hoomd.md.integrate.mode_standard(dt=dt);  
    all = hoomd.group.all();
    integrator = hoomd.md.integrate.nve(group=all);
    integrator.randomize_velocities(kT=0.8, seed=random.randint(1,10**4))

    ## Write output trajectory file
    hoomd.dump.gsd(trajectory_filename, period=steps/traj_frames, group=all, overwrite=True);


    ## Run the simulation
    hoomd.run(steps);

    ## Extract monomer coordinates
    traj =  gsd.hoomd.open(trajectory_filename)
    pos = traj[-1].particles.position

    return pos

def radius_gyration(pos):
    N = len(pos)
    x = []
    y = []
    z = []
    for k in np.arange(N):
        x.append(pos[k][0])
        y.append(pos[k][1])
        z.append(pos[k][2])

    r_mean = [np.mean(x),np.mean(y),np.mean(z)]
    rg_sq = 0
    for k in np.arange(N):
        rg_sq = rg_sq + sum((r_mean - pos[k])**2)

    rg_sq = rg_sq / N
    rg = np.sqrt(rg_sq)

    return rg

def render_polymer(pos, sigma = 1.0, color = [0.41, 0.30, 0.59]):
    scene = fresnel.Scene()
    N = len(pos)
    
    geometry = fresnel.geometry.Sphere(scene, N=N, position = pos, radius = sigma)
    geometry.material = fresnel.material.Material(color=fresnel.color.linear(color),
                                                roughness=0.8)

    scene.camera = fresnel.camera.Orthographic.fit(scene, view='isometric', margin=0.0)


    out = fresnel.preview(scene)
    image = PIL.Image.fromarray(out[:], mode='RGBA')
    image.save('polymer_render.png')

    return IPython.display.Image('polymer_render.png')

def randwalk2D(steps):
    X = []
    Y = []
    R = []
    x = 0
    y = 0
    for i in range(steps):
        X.append(x)
        Y.append(y)
        r = np.random.randint(1,5)
        R.append(r)
        if r == 1:
            x = x + 1
        elif r == 2:
            x = x - 1
        elif r == 3:
            y = y + 1
        elif r == 4:
            y = y - 1
    
    return X, Y

def Visualize2DRandomWalk(steps, Rg = 'on', end2end = 'on'):
    
    if steps >= 10000:
        linewidth = 1
        scale = 2
    else:
        linewidth = 2
        scale = 3
    X, Y = randwalk2D(steps)

    r_mean, rg = radius_gyration2DRW(X,Y,np.zeros(len(X)))

    fig, ax = plt.subplots(figsize=(10,10))

    xmin = r_mean[0] - scale*rg
    xmax = r_mean[0] + scale*rg
    ymin = r_mean[1] - scale*rg
    ymax = r_mean[1] + scale*rg

    plt.axis([xmin, xmax, ymin, ymax])
    plt.axis('Off')
    plt.plot(X,Y, 'k-', linewidth = linewidth)
    
    if Rg == 'on':
        color = np.array([112,48,160])/255
        draw_circle = mpl.patches.Circle((r_mean[0],r_mean[1]),rg, fill=True,edgecolor = color, facecolor = color, alpha = 0.5)
        ax.add_artist(draw_circle)
    if end2end == 'on':
        color = np.array([203,112,18])/255
        plt.plot([X[0], X[-1]],[Y[0],Y[-1]],'o--', color=color, linewidth = linewidth)
        
 def radius_gyration_2DRW(X,Y,Z):
    
    x = np.array(X)
    y = np.array(Y)
    z = np.zeros(len(X))
    N = len(x)
    pos = [x,y,z]
    pos = np.array(pos).transpose()

    r_mean = [np.mean(x),np.mean(y),np.mean(z)]
    rg_sq = 0
    
    for k in np.arange(N):
        rg_sq = rg_sq + sum((r_mean - pos[k])**2)

    rg_sq = rg_sq / N
    rg = np.sqrt(rg_sq)

    return r_mean, rg

        
