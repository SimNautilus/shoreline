void computeScaledJacobian(apf::Mesh2 * msh)
{
  // setup an apf::Field to apply displacements to
  apf::Field* j = m_msh->findField("scaledJacobian");
  if(!j) j = apf::createFieldOn(m_msh, "scaledJacobian", apf::DOUBLE);

  // get the global numbering of the mesh
  apf::GlobalNumbering* gn = m_msh->findGlobalNumbering("Global");
  // loop over vertices and apply displacement to the Field object
  apf::MeshIterator* it = m_msh->begin(0);
  apf::MeshEntity* v;
  while( (v = m_msh->iterate(it)) ) {
    long g = apf::getNumber(gn, v, 0, 0);
    double dx = dispArrays[0][g];
    double dy = dispArrays[1][g];
    double dz = dispArrays[2][g];
    apf::Vector3 disp(dx, dy, dz);
    setVector(d, v, 0, disp);
  }
  m_msh->end(it); // free the iterator
  apf::displaceMesh(m_msh, d);
  apf::synchronizeFieldData<double>(m_msh->getCoordinateField()->getData(),NULL);
  m_msh->acceptChanges();
}
