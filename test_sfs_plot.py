"""Generate a simulated SFS heatmap to demonstrate the updated plotting code."""
import tempfile
import numpy as np
from pathlib import Path


def make_folded_1d_sfs(n, total_snps=100000):
    sfs = np.zeros(n + 1)
    for i in range(1, n // 2 + 1):
        sfs[i] = total_snps / i
    sfs[0] = total_snps * 5
    sfs[1:n // 2 + 1] += np.random.poisson(500, size=n // 2)
    return sfs.astype(float)


def make_folded_2d_sfs(n1, n2, total_snps=50000, fst=0.1):
    sfs = np.zeros((n1 + 1, n2 + 1))
    for i in range(n1 + 1):
        for j in range(n2 + 1):
            if i + j == 0:
                sfs[i, j] = total_snps * 3
            elif i + j <= (n1 + n2) // 2:
                dist = np.sqrt(i ** 2 + j ** 2)
                if dist > 0:
                    sfs[i, j] = max(0, total_snps / (dist ** (1.5 + fst * 5)) + np.random.poisson(10))
    return sfs.astype(float)


def write_dadi_1d(path, sfs):
    n = len(sfs)
    with open(path, "w") as f:
        f.write(f"{n} unfolded\n")
        f.write(" ".join(f"{v:.1f}" for v in sfs) + "\n")
        mask = ["1"] + ["0"] * (n - 2) + ["1"]
        f.write(" ".join(mask) + "\n")


def write_dadi_2d(path, sfs):
    n1, n2 = sfs.shape
    with open(path, "w") as f:
        f.write(f"{n1} {n2} unfolded\n")
        vals = sfs.flatten()
        f.write(" ".join(f"{v:.1f}" for v in vals) + "\n")
        mask = np.zeros_like(sfs, dtype=int)
        mask[0, 0] = 1
        mask[-1, -1] = 1
        f.write(" ".join(str(int(v)) for v in mask.flatten()) + "\n")


def main():
    np.random.seed(42)
    pop_ids = ["P.trichocarpa", "P.balsamifera", "P.deltoides"]
    proj = 8

    with tempfile.TemporaryDirectory() as tmpdir:
        sfs_dir = Path(tmpdir)
        dadi_dir = sfs_dir / "dadi"
        dadi_dir.mkdir()

        for pop in pop_ids:
            sfs_1d = make_folded_1d_sfs(proj)
            write_dadi_1d(dadi_dir / f"{pop}-{proj}.sfs", sfs_1d)

        fst_map = {
            ("P.trichocarpa", "P.balsamifera"): 0.05,
            ("P.trichocarpa", "P.deltoides"): 0.20,
            ("P.balsamifera", "P.deltoides"): 0.15,
        }
        for i, pop1 in enumerate(pop_ids):
            for j, pop2 in enumerate(pop_ids):
                if i < j:
                    fst = fst_map.get((pop1, pop2), 0.1)
                    sfs_2d = make_folded_2d_sfs(proj, proj, fst=fst)
                    write_dadi_2d(dadi_dir / f"{pop1}-{pop2}.sfs", sfs_2d)

        from eprobe.popgen.assess import generate_multipop_sfs_heatmap

        output_path = Path("/Users/zh384/Desktop/scripts_dev/test_sfs_demo.png")
        result = generate_multipop_sfs_heatmap(
            sfs_dir=sfs_dir,
            pop_ids=pop_ids,
            output_path=output_path,
            title="Simulated SFS Demo (projection=8)",
        )
        if result.is_ok():
            print(f"Plot saved to: {result.unwrap()}")
        else:
            print(f"Error: {result.unwrap_err()}")


if __name__ == "__main__":
    main()
