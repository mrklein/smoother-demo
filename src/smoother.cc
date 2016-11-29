/*
 * smoother: small demonstration on how to:
 * - set initial field values
 * - smooth them
 * - reset them
 */

#include "argList.H"
#include "Hash.H"
#include "Time.H"
#include "volFields.H"
#include "fvc.H"


int
main(int argc, char **argv)
{
  using namespace Foam;

  argList::noBanner();
  argList args{argc, argv};

  if (!args.checkRootCase()) {
    Foam::FatalError.exit();
  }

  Info << "PID: " << pid() << nl << endl;

  Time run_time{Time::controlDictName, args};
  fvMesh mesh({fvMesh::defaultRegion,
              run_time.timeName(),
              run_time,
              IOobject::MUST_READ});

  volScalarField alpha({"alpha", run_time.timeName(), mesh, IOobject::NO_READ},
                       mesh,
                       {"0", dimless, 0.0});

  // Generating initial values
  Info << "#" << run_time.timeIndex() << endl;
  forAll(alpha, i) {
    const auto& Ci = mesh.C()[i];
    const auto r = ::sqrt(sqr(Ci.x()) + sqr(Ci.y()));
    const auto s = 2*::j1(M_PI*r)/M_PI/r;
    if (Ci.z() < s)
      alpha.ref()[i] = 1;
  }

  alpha.write();

  // Smoothing
  run_time++;
  Info << "#" << run_time.timeIndex() << endl;
  const auto nSmoothIter = 2;
  for (auto i = 0; i < nSmoothIter; i++) {
    surfaceScalarField alphaf{fvc::interpolate(alpha)};
    alpha = fvc::average(alphaf);
    // Or just
    // alpha = fvc::average(fvc::interpolate(alpha));
  }

  alpha.write();

  // Resetting field to 0
  run_time++;
  Info << "#" << run_time.timeIndex() << endl;
  alpha.ref() = 0.0;
  alpha.write();

  return 0;
}

