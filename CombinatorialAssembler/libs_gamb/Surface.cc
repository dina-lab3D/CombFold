#include "Shou.h"
#include "Surface.h"
#include "Atom.h"

Surface::Surface() : Molecule<SurfacePoint>()
{}

int Surface::readShouFile(std::istream& shouStream) {
  return readSurfaceFile(shouStream);
}

int Surface::readSurfaceFile(std::istream& inFile) {
  if(isBinaryFile(inFile)) {
    //cerr << "BINARY" << endl;
    return readBinaryFile(inFile);
  } else {
    //cerr << "ASCII" << endl;
    return readASCIIFile(inFile);
  }
}


int Surface::readSurfaceFile(const std::string fileName) {
  std::ifstream inFile(fileName.c_str());
  if (!inFile) {
    std::cerr << "Can't open surface file: " << fileName << std::endl;
    exit(1);
  }
  return readSurfaceFile(inFile);
}

void Surface::output2PDB(const std::string outFileName) {
  std::ofstream outFile(outFileName.c_str());
  unsigned int counter = 0;
  for(iterator it = begin(); it != end() ; it++, counter++) {
    SurfacePoint::SurfaceType nodeType = (*it).surfaceType();
     unsigned int resId = 10;
     switch(nodeType) {
      case SurfacePoint::Caps:  resId = 10; break;
      case SurfacePoint::Pits:  resId = 20; break;
      case SurfacePoint::Belts: resId = 30; break;
      default: break;
     }
     Atom atom(it->position(), ' ', counter, resId, "PRB", 'X');
     outFile << atom << std::endl;
  }
}

int Surface::readBinaryFile(std::istream& inFile) {
  int atoms[3];
  float coord[3];
  float area;
  float normal[3];

  while(!inFile.eof()) {
    // read 3 atom numbers
    for(int i=0; i<3; i++) inFile.read((char *)&atoms[i], sizeof(int));

    // read point coordinates
    for(int i=0; i<3; i++) inFile.read((char *)&coord[i], sizeof(float));

    // read area
    inFile.read((char *)&area, sizeof(float));

    // read normal
    for(int i=0; i<3; i++) inFile.read((char *)&normal[i], sizeof(float));

    if(inFile.eof()) break;

    Vector3 point(coord);
    Vector3 norm(normal);
    SurfacePoint sp(point, norm, area, atoms[0]-1, atoms[1]-1, atoms[2]-1);
    push_back(sp);
  }
  std::cerr << "Surface size: " << size() << std::endl;
  return size();
}

int Surface::readASCIIFile(std::istream& shouStream) {
  while (!shouStream.eof())  {
    std::string line;
    getline(shouStream, line);
    if (Shou::isRec(line.c_str()))
      add(SurfacePoint(line.c_str()));
  }
  return size();
}

float Surface::area() const {
  float area=0;
  for (const_iterator it=begin(); it != end(); ++it) {
    area+= it->surfaceArea();
  }
  return area;
}

void Surface::outputVRML(std::ostream& VRMLstream,
                         const float red,
                         const float green,
                         const float blue) const
{
  const char header [] = "#VRML V1.0 ascii";  // VRML 1.0 file header
  const char separator[] = "Separator {\n";   // VRML spearator
  const char ind[] = "  ";                    // indentation

  VRMLstream << header << '\n' << '\n';
  VRMLstream << separator;
  VRMLstream << ind << "Material {\n";
  VRMLstream << ind << ind << "diffuseColor "
             << red << ' ' << green << ' ' << blue << '\n';
  VRMLstream << ind << "}\n";

  VRMLstream << ind << "Coordinate3 {\n";
  VRMLstream << ind << ind << "point [\n";
  for (const_iterator it=begin(); it != end(); ++it) {
    if (it != begin())
      VRMLstream << ',' << '\n';
    const Vector3 pos = (*it).position();
    const Vector3 tip = pos + (*it).normal();
    VRMLstream << ind << ind << ind
               << pos[0] << ' ' << pos[1] << ' ' << pos[2] << ',' << '\n';
    VRMLstream << ind << ind << ind
               << tip[0] << ' ' << tip[1] << ' ' << tip[2];
  }
  VRMLstream << '\n';
  VRMLstream << ind << ind << ']' << '\n';
  VRMLstream << ind << '}' << '\n';
  VRMLstream << ind << "IndexedLineSet {\n";
  VRMLstream << ind << ind << "coordIndex [\n";

  for (unsigned i=0; i<size(); ++i) {
    VRMLstream << ind << ind << ind << i*2 << ", " << i*2+1 << ", " << -1;
    if (i != size())
      VRMLstream << ',';
    VRMLstream << '\n';
  }

  VRMLstream << ind << ind << ']' << '\n';
  VRMLstream << ind << '}' << '\n';
  VRMLstream << '}' << '\n';
}

std::ostream& operator<<(std::ostream& s, const Surface& surf)
{
  for (Surface::const_iterator it=surf.begin(); it != surf.end() ; ++it)
    s << *it ;
  s << std::endl;
  return s;
}

bool Surface::isBinaryFile(std::istream& inFile) {
  int TEST_NUM = 100;
  char ch;
  for(int i=0; i<TEST_NUM; i++) {
    inFile.get(ch);
    if(!isascii(ch)) {
      inFile.seekg(0, std::ios::beg);
      return true;
    }
  }

  inFile.seekg(0, std::ios::beg);
  return false;
}
