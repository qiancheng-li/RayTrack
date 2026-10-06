import os
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry
from pyproj import Transformer
from scipy.io import loadmat
from scipy.ndimage import gaussian_filter

from utils import *

def visualize_candidate_region(valid_xy, poss_valid_idx, path_xy, actual_xy, start_xy=None, end_xy=None, title="Candidate Region"):
    valid_xy = np.asarray(valid_xy)
    poss_valid_idx = np.asarray(poss_valid_idx, dtype=int)

    selected_xy = valid_xy[poss_valid_idx]

    plt.figure(figsize=(7, 6))

    # All valid points.
    plt.scatter(valid_xy[:, 0], valid_xy[:, 1], s=8, alpha=0.25, label="All valid points")

    # Selected candidate points.
    plt.scatter(selected_xy[:, 0], selected_xy[:, 1], s=12, alpha=0.9, label="Selected candidates")

    # Start and end points.
    if start_xy is not None:
        plt.scatter(start_xy[0], start_xy[1], s=80, marker="o", label="Start")
    if end_xy is not None:
        plt.scatter(end_xy[0], end_xy[1], s=80, marker="x", label="End")

    plt.scatter(path_xy[:, 0], path_xy[:, 1], s=10, alpha=1, label="path")
    plt.scatter(actual_xy[:, 0], actual_xy[:, 1], s=10, alpha=1, label="actual")

    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title(title)
    plt.axis("equal")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

# =========================================================
# basic io
# =========================================================
def _load_mat_variable(mat_path):
    data = loadmat(mat_path)
    for k, v in data.items():
        if not k.startswith("__"):
            return v
    raise ValueError(f"No valid variable found in {mat_path}")


# =========================================================
# coordinate transform
# NOTE:
# Convert latitude/longitude to the local metric coordinate system used by the radio map.
# =========================================================
def sionna_coord_batch(latlon_array):
    """
    latlon_array: [N,2] -> [[lat, lon], ...]
    return: [N,2] -> [[x, y], ...]
    Return local metric coordinates centered on the scene bounding box.
    """
    latlon_array = np.asarray(latlon_array, dtype=float)
    if latlon_array.ndim == 1:
        latlon_array = latlon_array.reshape(1, 2)

    assert latlon_array.shape[1] == 2, "Input must be Nx2 [lat, lon]"

    bl = (-122.312299126737, 47.65074792715546)
    ur = (-122.30187970485449, 47.65690207846001)
    box = [bl, (bl[0], ur[1]), ur, (ur[0], bl[1])]

    epsg_code = get_utm_epsg_code_from_gps(bl[0], bl[1])
    to_utm = Transformer.from_crs("EPSG:4326", epsg_code, always_xy=True)

    ground_utm = np.array([to_utm.transform(lon, lat) for lon, lat in box])
    poly = shapely.geometry.Polygon(ground_utm)
    center_x = poly.envelope.centroid.x
    center_y = poly.envelope.centroid.y

    lats = latlon_array[:, 0]
    lons = latlon_array[:, 1]
    utm_x, utm_y = to_utm.transform(lons, lats)

    x = np.asarray(utm_x) - center_x
    y = np.asarray(utm_y) - center_y
    return np.stack([x, y], axis=1)


def xy_distance(x1, y1, x2, y2):
    return np.hypot(np.asarray(x1) - np.asarray(x2), np.asarray(y1) - np.asarray(y2))


# =========================================================
# xy -> powermat row/col
# Direct XY-to-index conversion.
# =========================================================
def xy_to_rc(xy, powermat_shape):
    """
    Directly interpret xy as powermat index coordinates.

    Parameters
    ----------
    xy : [2] or [N,2]
    powermat_shape : (H, W)

    Returns
    -------
    rc : [2] or [N,2], int
        row/col after rounding and clipping
    """
    xy = np.asarray(xy, dtype=float)
    H, W = powermat_shape

    single = False
    if xy.ndim == 1:
        xy = xy.reshape(1, 2)
        single = True

    # Legacy convention: x maps to rows and y maps to columns.
    row = np.round(xy[:, 0]).astype(int)
    col = np.round(xy[:, 1]).astype(int)

    row = np.clip(row, 0, H - 1)
    col = np.clip(col, 0, W - 1)

    rc = np.stack([row, col], axis=1)
    return rc[0] if single else rc


# =========================================================
# candidate region
# =========================================================
def find_possible_xy(start_xy, end_xy, all_valid_xy, expand_factor=2.0):
    start_xy = np.asarray(start_xy, dtype=float).reshape(2,)
    end_xy = np.asarray(end_xy, dtype=float).reshape(2,)
    all_valid_xy = np.asarray(all_valid_xy, dtype=float)

    dis1 = np.linalg.norm(all_valid_xy - start_xy[None, :], axis=1)
    dis2 = np.linalg.norm(all_valid_xy - end_xy[None, :], axis=1)
    total_dis = dis1 + dis2
    pnt_dis = np.linalg.norm(start_xy - end_xy)

    return np.where(total_dis < expand_factor * pnt_dis)[0]

def xy_to_rc_from_cell_centers(xy, cell_centers, map_shape):
    """
    Use cell_centers to convert XY coordinates to nearest grid row/col
    by axis rounding (no KD-tree).

    Parameters
    ----------
    xy : [2] or [N,2]
    cell_centers : [H,W,2 or 3]
    map_shape : (H,W)

    Returns
    -------
    rc : [2] or [N,2], int
    """
    xy = np.asarray(xy, dtype=float)
    single = False
    if xy.ndim == 1:
        xy = xy.reshape(1, 2)
        single = True

    H, W = map_shape
    # print(cell_centers.shape)

    # Extract the x axis from columns and the y axis from rows.
    x_axis = cell_centers[0, :, 0]
    y_axis = cell_centers[:, 0, 1]

    dx = np.mean(np.diff(x_axis))
    dy = np.mean(np.diff(y_axis))

    x0 = x_axis[0]
    y0 = y_axis[0]

    row = np.round((xy[:, 0] - x0) / dx).astype(int)
    col = np.round((xy[:, 1] - y0) / dy).astype(int)

    row = np.clip(row, 0, H - 1)
    col = np.clip(col, 0, W - 1)

    rc = np.stack([row, col], axis=1)
    return rc[0] if single else rc


# =========================================================
# sampling
# Sample the power map at the nearest grid cells.
# =========================================================
def find_power_from_xy_direct(path_xy, cell_centers, powermat):
    """
    Direct round-based sampling using cell_centers axes.
    """
    path_xy = np.asarray(path_xy, dtype=float)
    powermat = np.asarray(powermat, dtype=float)

    rc = xy_to_rc_from_cell_centers(path_xy, cell_centers, powermat.shape)
    rows = rc[:, 0]
    cols = rc[:, 1]
    return powermat[rows, cols]

def plot_find_path_stage(
    cell_centers,
    powermat,
    peak_mat,
    start_xy,
    end_xy,
    stage_path_xy=None,
    valid_xy=None,
    poss_valid_idx=None,
    title="find_path stage",
    show_powermat=False,
):
    """
    Visualize one stage inside find_path4_xy_sparse.

    Parameters
    ----------
    cell_centers : [H,W,2 or 3]
        Grid cell centers in XY coordinates.
    powermat : [H,W]
        Map values. Optional as background.
    peak_mat : [H,W]
        Peak map after 1/mat and gaussian smoothing.
    start_xy, end_xy : [2]
        Start and end points in XY.
    stage_path_xy : [L,2] or None
        Final path for this stage (e.g. the path returned at depth==1,
        or the candidate path corresponding to the chosen midpoint).
    valid_xy : [N,2] or None
        All valid points.
    poss_valid_idx : [K] or None
        Candidate subset indices into valid_xy.
    title : str
        Plot title.
    show_powermat : bool
        Whether to also show powermat as a faint background.
    """
    cell_centers = np.asarray(cell_centers, dtype=float)
    peak_mat = np.asarray(peak_mat, dtype=float)
    powermat = np.asarray(powermat, dtype=float)
    start_xy = np.asarray(start_xy, dtype=float).reshape(2,)
    end_xy = np.asarray(end_xy, dtype=float).reshape(2,)

    H, W = peak_mat.shape

    # Grid coordinates.
    X = cell_centers[:, :, 0]
    Y = cell_centers[:, :, 1]

    plt.figure(figsize=(16, 14))

    # Optional power-map background.
    if show_powermat:
        plt.pcolormesh(X, Y, powermat, shading="auto", alpha=0.25)

   

    # All valid points.
    if valid_xy is not None:
        valid_xy = np.asarray(valid_xy, dtype=float)
        plt.scatter(
            valid_xy[:, 0], valid_xy[:, 1],
            s=6, alpha=0.15, label="all valid"
        )

    # Candidate points.
    if valid_xy is not None and poss_valid_idx is not None and len(poss_valid_idx) > 0:
        poss_valid_idx = np.asarray(poss_valid_idx, dtype=int)
        cand_xy = valid_xy[poss_valid_idx]
        plt.scatter(
            cand_xy[:, 0], cand_xy[:, 1],
            s=10, alpha=0.7, label="candidate"
        )

    # Peak-map background.
    pcm = plt.pcolormesh(X, Y, peak_mat, shading="auto", alpha=0.85)
    plt.colorbar(pcm, label="peak_mat")

    # Start and end points.
    plt.scatter(start_xy[0], start_xy[1], s=90, marker="o", label="start")
    plt.scatter(end_xy[0], end_xy[1], s=90, marker="x", label="end")

    # Path selected at the current stage.
    if stage_path_xy is not None:
        stage_path_xy = np.asarray(stage_path_xy, dtype=float)
        plt.plot(
            stage_path_xy[:, 0], stage_path_xy[:, 1],
            "-r", linewidth=2, label="stage path"
        )
        plt.scatter(
            stage_path_xy[:, 0], stage_path_xy[:, 1],
            s=14, c="r", alpha=0.8
        )

    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title(title)
    plt.axis("equal")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.show()


# =========================================================
# recursive path finder
# Candidate points use the same XY coordinate system as cell_centers.
# =========================================================
def find_path4_xy_sparse(
    poss_valid_idx,   # indices into valid_xy
    start_xy,
    end_xy,
    exp_data,
    depth,
    valid_xy,
    powermat,
    cell_centers,
    smooth_sigma=6,
):
    """
    Estimate a path by evaluating candidate midpoints on the radio-map grid.
    """
    poss_valid_idx = np.asarray(poss_valid_idx, dtype=int).reshape(-1)
    start_xy = np.asarray(start_xy, dtype=float).reshape(2,)
    end_xy = np.asarray(end_xy, dtype=float).reshape(2,)
    exp_data = np.asarray(exp_data, dtype=float).reshape(-1)
    valid_xy = np.asarray(valid_xy, dtype=float)
    powermat = np.asarray(powermat, dtype=float)
    cell_centers = np.asarray(cell_centers, dtype=float)

    desired_len = len(exp_data)

    if desired_len == 1:
        return start_xy.reshape(1, 2)

    if desired_len == 2:
        return np.vstack([start_xy, end_xy])

    H, W = powermat.shape
    mat = np.full((H, W), np.inf, dtype=float)

    for vidx in poss_valid_idx:
        mid_xy = valid_xy[vidx]

        # Map the candidate midpoint to its nearest grid cell.
        r, c = xy_to_rc_from_cell_centers(mid_xy, cell_centers, powermat.shape)

        start_distance = np.linalg.norm(start_xy - mid_xy)
        end_distance = np.linalg.norm(end_xy - mid_xy)

        denom = start_distance + end_distance
        if denom == 0:
            first_len = desired_len // 2
        else:
            first_len = int(np.floor(desired_len * start_distance / denom))
        second_len = desired_len - first_len + 1

        standard_speed = 0.2624
        speed = (start_distance + end_distance) / max(desired_len - 1, 1)
        if speed == 0:
            penalty_factor = np.inf
        else:
            relative_speed = speed / standard_speed
            penalty_factor = max(relative_speed, 1.0 / relative_speed)

        x1 = np.linspace(start_xy[0], mid_xy[0], first_len)
        x2 = np.linspace(mid_xy[0], end_xy[0], second_len)
        y1 = np.linspace(start_xy[1], mid_xy[1], first_len)
        y2 = np.linspace(mid_xy[1], end_xy[1], second_len)

        xs = np.concatenate([x1, x2[1:]])
        ys = np.concatenate([y1, y2[1:]])
        path_xy = np.stack([xs, ys], axis=1)

        if len(path_xy) != desired_len:
            continue

        ss = find_power_from_xy_direct(path_xy, cell_centers, powermat)
        candidate_loss = np.sum(np.abs(ss - exp_data)) * penalty_factor
        # plt.figure()
        # plt.plot(exp_data)
        # plt.plot(ss)
        # plt.show()

        # Store the candidate loss at the corresponding grid cell.
        mat[c, r] = candidate_loss

    peak_mat = 1.0 / mat
    peak_mat[~np.isfinite(peak_mat)] = 0.0
    peak_mat = gaussian_filter(peak_mat, sigma=smooth_sigma)

    

    best_flat = int(np.argmax(peak_mat))
    best_r, best_c = np.unravel_index(best_flat, peak_mat.shape)

    # Recover the physical XY coordinate of the best grid cell.
    mid_xy = np.array(cell_centers[best_r, best_c, :2], dtype=float)
    # print(start_xy, end_xy, mid_xy)
    # plt.figure()
    # plt.imshow(peak_mat)
    # plt.title("peak_mat")
    # plt.colorbar()
    # plt.show()

    start_distance = np.linalg.norm(start_xy - mid_xy)
    end_distance = np.linalg.norm(end_xy - mid_xy)
    denom = start_distance + end_distance
    if denom == 0:
        first_len = desired_len // 2
    else:
        first_len = int(np.floor(desired_len * start_distance / denom))
    second_len = desired_len - first_len + 1

    if depth == 1:
        x1 = np.linspace(start_xy[0], mid_xy[0], first_len)
        x2 = np.linspace(mid_xy[0], end_xy[0], second_len)
        y1 = np.linspace(start_xy[1], mid_xy[1], first_len)
        y2 = np.linspace(mid_xy[1], end_xy[1], second_len)

        xs = np.concatenate([x1, x2[1:]])
        ys = np.concatenate([y1, y2[1:]])

        stage_path_xy = np.stack([xs, ys], axis=1)

        # plot_find_path_stage(
        #     cell_centers=cell_centers,
        #     powermat=powermat,
        #     peak_mat=peak_mat,
        #     start_xy=start_xy,
        #     end_xy=end_xy,
        #     stage_path_xy=stage_path_xy,
        #     valid_xy=valid_xy,
        #     poss_valid_idx=poss_valid_idx,
        #     title=f"find_path stage | depth={depth} | len={desired_len}",
        #     show_powermat=False,
        # )
        return np.stack([xs, ys], axis=1)

    path1 = find_path4_xy_sparse(
        poss_valid_idx,
        start_xy,
        mid_xy,
        exp_data[:first_len],
        depth - 1,
        valid_xy,
        powermat,
        cell_centers,
        smooth_sigma=smooth_sigma,
    )

    path2 = find_path4_xy_sparse(
        poss_valid_idx,
        mid_xy,
        end_xy,
        exp_data[first_len - 1:],
        depth - 1,
        valid_xy,
        powermat,
        cell_centers,
        smooth_sigma=smooth_sigma,
    )

    return np.vstack([path1[:-1], path2])

def load_height_map(height_tag, base_save_dir):
    """Load the tracking radio map for a transmitter-height setting."""
    npz_path = Path(base_save_dir) / f"lab1_tracking_radio_map_{height_tag}.npz"
    data = np.load(npz_path)

    powermat = data["rss_db"]
    cell_centers = data["cell_centers"]
    return powermat, cell_centers

def load_LoS_map(height_tag, base_save_dir="./height_sensitivity_results"):
    npz_path = os.path.join(base_save_dir, f"{height_tag}.npz")
    data = np.load(npz_path)

    powermat = data["rss_db"]
    cell_centers = data["cell_centers"]
    return powermat, cell_centers 

# =========================================================
# top-level tester
# =========================================================
def algorithm_tester_xy_sparse(datafilename, samplespersearch, depth, height_tag):
    repo_root = Path(__file__).resolve().parents[1]
    powergridname = repo_root / "data" / "algorithm_inputs" / "lab1_radio_map_points.mat"

    slope = 0.518
    offset = -13.243

    experimental_data = np.loadtxt(datafilename, delimiter=",", skiprows=1)
    experimental_data = np.asarray(experimental_data, dtype=float)

    rssi_data = experimental_data[:, 2] * slope + offset

    power_map_grid = _load_mat_variable(powergridname)
    print(power_map_grid.shape)
    print(power_map_grid[0])

    radio_map_dir = repo_root / "data" / "algorithm_inputs"
    powermat, cell_centers = load_height_map(height_tag, base_save_dir=radio_map_dir)
    powermat = np.asarray(powermat, dtype=float)
    cell_centers = np.asarray(cell_centers, dtype=float)

    power_map_grid = np.asarray(power_map_grid, dtype=float)
    powermat = np.asarray(powermat, dtype=float)

    # Retain valid XY candidates; no separate row/column array is required.
    valid_latlon = power_map_grid[:, :2]
    valid_xy = sionna_coord_batch(valid_latlon)

    actual_xy = sionna_coord_batch(experimental_data[:, :2])

    final_path_xy = np.zeros((len(rssi_data), 2), dtype=float)

    chunkNum = int(np.floor(len(rssi_data) / samplespersearch))
    time_mat = np.zeros((chunkNum, 1), dtype=float)

    for chunkID in range(1, chunkNum + 1):
        print(chunkID, chunkNum)

        start_row = (chunkID - 1) * samplespersearch
        start_xy = actual_xy[start_row]

        if chunkID * samplespersearch < actual_xy.shape[0]:
            end_xy = actual_xy[chunkID * samplespersearch]
        else:
            end_xy = actual_xy[chunkID * samplespersearch - 1]

        poss_valid_idx = find_possible_xy(start_xy, end_xy, valid_xy)
        # visualize_candidate_region(
        #     valid_xy,
        #     poss_valid_idx,
        #     start_xy=start_xy,
        #     end_xy=end_xy,
        #     title=f"Chunk {chunkID} Candidate Region"
        # )

        t0 = time.perf_counter()

        if chunkID * samplespersearch < actual_xy.shape[0]:
            rss_chunk = rssi_data[
                (chunkID - 1) * samplespersearch : chunkID * samplespersearch + 1
            ]

            path_xy = find_path4_xy_sparse(
                poss_valid_idx,
                start_xy,
                end_xy,
                rss_chunk,
                depth,
                valid_xy,
                powermat,
                cell_centers
            )
            final_path_xy[
                (chunkID - 1) * samplespersearch : chunkID * samplespersearch, :
            ] = path_xy[:-1]

        else:
            idx_start = (chunkID - 1) * samplespersearch
            idx_end = chunkID * samplespersearch

            rss_chunk = rssi_data[
                (chunkID - 1) * samplespersearch : chunkID * samplespersearch
            ]

            path_xy = find_path4_xy_sparse(
                poss_valid_idx,
                start_xy,
                end_xy,
                rss_chunk,
                depth,
                valid_xy,
                powermat,
                cell_centers
            )
            final_path_xy[idx_start:idx_end, :] = path_xy

        time_mat[chunkID - 1, 0] = time.perf_counter() - t0

    if chunkNum > 0 and chunkID * samplespersearch < actual_xy.shape[0]:
        start_xy = actual_xy[chunkID * samplespersearch]
        end_xy = actual_xy[-1]
        poss_valid_idx = find_possible_xy(start_xy, end_xy, valid_xy)

        path_xy = find_path4_xy_sparse(
            poss_valid_idx,
            start_xy,
            end_xy,
            rssi_data[chunkID * samplespersearch :],
            depth,
            valid_xy,
            powermat,
            cell_centers
        )
        final_path_xy[chunkID * samplespersearch :, :] = path_xy

    elif chunkNum == 0:
        start_xy = actual_xy[0]
        end_xy = actual_xy[-1]
        poss_valid_idx = find_possible_xy(start_xy, end_xy, valid_xy)

        path_xy = find_path4_xy_sparse(
            poss_valid_idx,
            start_xy,
            end_xy,
            rssi_data,
            depth,
            valid_xy,
            powermat,
            cell_centers
        )
        final_path_xy[:, :] = path_xy

    # visualize_candidate_region(
    #         valid_xy,
    #         poss_valid_idx,
    #         final_path_xy,
    #         actual_xy,
    #         start_xy=start_xy,
    #         end_xy=end_xy,
    #         title=f"Chunk {chunkID} Candidate Region"
    #     )

    ERR_MAT = np.linalg.norm(final_path_xy - actual_xy, axis=1)

    errors_np = np.array(ERR_MAT, dtype=float)

    # 95% confidence interval of mean MAE
    n = len(errors_np)
    if n > 1:
        ci95 = 1.96 * np.std(errors_np, ddof=1) / np.sqrt(n)
    else:
        ci95 = 0.0

    mean_err = np.median(ERR_MAT)
    std_err = np.std(ERR_MAT)

    return mean_err, std_err, ci95, ERR_MAT, time_mat, final_path_xy, actual_xy


# =========================================================
# main
# =========================================================
def main():
    repo_root = Path(__file__).resolve().parents[1]
    trajectory_dir = (
        repo_root
        / "data"
        / "experimental_standardized"
        / "Station-Lab1"
        / "lab1 data"
    )
    trajectory_names = [
        "11-14-1-interpolated.csv",
        "11-14-2-interpolated.csv",
        "11-14-3-interpolated.csv",
        "11-14-4-interpolated.csv",
        "11-14-5-interpolated.csv",
        "11-14-6-interpolated.csv",
        "11-14-7-interpolated.csv",
        "11-26-8-interpolated.csv",
        "11-26-9-interpolated.csv",
        "11-26-10-interpolated.csv",
    ]
    datafiles = [trajectory_dir / name for name in trajectory_names]

    samplespersearch = 90
    depth = 1

    height_tag = "13m"

    mean_errs = []
    std_errs = []
    mae_ci_results = []
    all_errs = []

    for datafile in datafiles:
        mean_err, std_err, ci95, ERR_MAT, time_mat, final_path_xy, actual_xy = algorithm_tester_xy_sparse(
            datafilename=datafile,
            samplespersearch=samplespersearch,
            depth=depth,
            height_tag=height_tag,
        )
        mean_errs.append(mean_err)
        std_errs.append(std_err)
        mae_ci_results.append(ci95)
        all_errs.append(ERR_MAT)

        print("========== Result ==========")
        print(f"Height:          {height_tag}")
        print(f"Median error (m): {mean_err:.4f}")
        print(f"Std error (m):    {std_err:.4f}")
        if len(time_mat) > 0:
            print(f"Mean time/chunk:  {np.mean(time_mat):.4f} s")
        else:
            print("Mean time/chunk:  N/A")
        print(f"Num samples:      {len(ERR_MAT)}")
        print()

    print("========== Overall Result ==========")
    print(f"Median error (m): {np.median(np.concatenate(all_errs)):.4f}")
    print(f"Std error (m):    {np.std(np.concatenate(all_errs)):.4f}")
    print()


if __name__ == "__main__":
    main()
