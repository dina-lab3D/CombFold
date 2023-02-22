
#ifndef ENERGY_ATOM_H
#define ENERGY_ATOM_H

class EnergyAtom{
 public:

    virtual float getRadius() const = 0;
    virtual float getEpsilon() const = 0;
    virtual bool isHydrogen() const = 0;
    virtual float getCharge() const = 0;
    virtual bool isDonor() const = 0;
    virtual bool isAcceptor() const = 0;
    virtual const Vector3& getHBDirection() const = 0;
    virtual Vector3 position() const = 0;

    virtual ~EnergyAtom() {}
};

#endif

