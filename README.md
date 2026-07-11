# The GAP package CanonicalPcPres

CanonicalPcPres computes canonical polycyclic presentations of finite solvable
groups. For a finite solvable pc group `G`, it returns the canonical pc group
`can(G)` together with an explicit isomorphism from `G` to `can(G)`.

The package is a proof-of-concept implementation of the algorithm described in
the accompanying paper, [*Computing canonical labellings of finite solvable
groups*](https://arxiv.org/abs/2606.26030), by Santiago Barrera Acevedo,
Heiko Dietrich, and Max Horn.

## Installation

CanonicalPcPres requires GAP 4.13 or newer and the GAP packages
[ANUPQ](https://gap-packages.github.io/anupq/) and
[orb](https://gap-packages.github.io/orb/) and
[images](https://gap-packages.github.io/images/).

Install a release archive in a GAP package directory, or clone the repository
there:

```sh
cd path/to/gap/pkg
git clone https://github.com/gap-packages/CanonicalPcPres
```

Then load the package in GAP:

```gap
LoadPackage("CanonicalPcPres");
```

## Example

```gap
gap> G := SmallGroup(24*7, 54);;
gap> iso := IsomorphismCanonicalPcGroup(G);
[ f3, f2, f1, f4, f5 ] -> [ f1, f2, f3, f4, f5 ]
gap> Source(iso) = G;
true
gap> IsBijective(iso);
true
gap> Range(iso) = CanonicalPcGroup(G);
true
```

The complete manual is available on the
[package website](https://gap-packages.github.io/CanonicalPcPres/). After
loading the package, it is also available through GAP's help system with
`?CanonicalPcPres`.

## Contact

Questions and bug reports should be submitted through the
[issue tracker](https://github.com/gap-packages/CanonicalPcPres/issues).
The authors and maintainers are Santiago Barrera Acevedo, Heiko Dietrich, and
Max Horn; their contact details are listed in `PackageInfo.g`.

## License

CanonicalPcPres is distributed under the GNU General Public License,
version 2 or later. See [LICENSE](LICENSE).
