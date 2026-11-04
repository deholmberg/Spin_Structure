#include "Particle.h"
using namespace std;

//----------------------------------------------------------------------------------------Constructeurs--------------------------------------------------------------------//
PARTICLE::PARTICLE() : momentum(), pid()
{}
PARTICLE::PARTICLE(TVector3 momentum, TVector3 position, int pid, int chi2pid, int status, int charge) : momentum(momentum), position(position), pid(pid), chi2pid(chi2pid), status(status), charge(charge)
{}

PARTICLE::PARTICLE(TVector3 momentum, TVector3 position, int pid, int chi2pid, int status, int charge, double beta) : momentum(momentum), position(position), pid(pid), chi2pid(chi2pid), status(status), charge(charge), beta(beta)
{}


PARTICLE::PARTICLE(TVector3 momentum, int pid) : momentum(momentum), pid(pid)
{}

//----------------------------------------------------------------------------------------Setters----------------------------------------------------------------------------------------//
void PARTICLE::setmomentum(TVector3 momentum)
{
this->momentum = momentum;
}

void PARTICLE::setposition(TVector3 position)
{
  this->position = position;
}

void PARTICLE::setpid(int pid)
{
  this->pid = pid;
}

void PARTICLE::setchi2pid(int chi2pid)
{
  this->chi2pid = chi2pid;
}

void PARTICLE::setstatus(int status)
{
  this->status = status;
}

void PARTICLE::setcharge(int charge)
{
  this->charge = charge;
}

//----------------------------------------------------------------------------------------Getters----------------------------------------------------------------------------------------//
TVector3 PARTICLE::getmomentum() const
{
  return this->momentum;
}

TVector3 PARTICLE::getposition() const
{
  return this->position;
}

int PARTICLE::getpid() const
{
  return this->pid;
}

int PARTICLE::getchi2pid() const
{
  return this->chi2pid;
}

int PARTICLE::getstatus() const
{
  return this->status;
}

int PARTICLE::getcharge() const
{
  return this->charge;
}

//----------------------------------------------------------------------------------------Other methods----------------------------------------------------------------------------------------//

void PARTICLE::get_mass()
{
  if(this->pid==11) this->mass=0.51099895000/1000.;
  else if(this->pid==2212) this->mass=938.27208816/1000.;
  else cout << "PID unknown : " << this->pid << endl;
}

void PARTICLE::compute_energy()
{
  double E =  this->momentum.Mag2() + this->mass*this->mass;
  this->energy = TMath::Sqrt(E);
}

void PARTICLE::compute_theta()
{
  this->theta = (180. / TMath::Pi()) * this->momentum.Theta();
}

void PARTICLE::compute_phi()
{
  this->phi = (180. / TMath::Pi()) * this->momentum.Phi();
}

void PARTICLE::compute_LorentzVector()
{
  TLorentzVector V(this->momentum, this->energy);
  this->LorentzVector = V;
}

void PARTICLE::construct_all_properties()
{
  this->get_mass();
  this->compute_energy();
  this->compute_LorentzVector();
  this->compute_phi();
  this->compute_theta();
}
