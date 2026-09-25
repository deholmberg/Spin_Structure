#ifndef PARTICLE_MASSES_H
#define PARTICLE_MASSES_H

const double p_mass = 0.93827208816; // proton mass in GeV/c^2 from NIST
const double n_mass = 0.93956542052; // neutron mass in GeV/c^2 from NIST
const double nucleon_mass = (p_mass + n_mass)/2.0;
const double ppm_mass=0.13957039; // +-pion mass in GeV/c^2 from Particle Data Group
const double p0_mass= 4.5936/1000.0; // neutral pion mass in GeV/c^2 from Particle Data Group
const double k_mass = 0.493677; // K+- mass in GeV/c^2 from Particle Data Group
const double e_mass = 0.51099895/1000.0; // electron mass in GeV/c^2 from NIST
const double d_mass = 1.87561294257; // deuteron mass in GeV/c^2 from NIST
const double mu_p = 2.79; // Proton magnetic moment

const double beam_energy1 = 10.5473; // electron beam energy in GeV for RGC runs below 17065
const double beam_energy2 = 10.5563; // electron beam energy in GeV for RGC runs between 17067 to 17704
const double beam_energy3 = 10.5593; // electron beam energy in GeV for RGC runs above 17720

#endif
