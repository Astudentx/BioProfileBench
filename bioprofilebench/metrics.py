"""Benchmark metric calculations."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr, spearmanr


CORE_COLUMNS = [
    "Sample",
    "Taxonomic Level",
    "Info",
    "TP",
    "FP",
    "FN",
    "TN",
    "Precision",
    "Recall",
    "F1",
    "Spearman",
    "AbundancePrecision",
    "AbundanceRecall",
    "AbundanceF1",
    "L1",
    "L2",
    "BC",
    "rJSD",
]

EXTRA_COLUMNS = [
    "Jaccard",
    "WeightedJaccard",
    "Pearson",
    "Cosine",
    "AUPRC",
    "Accuracy",
    "Specificity",
    "NegativePredictiveValue",
    "MatthewsCorrelationCoefficient",
    "MeanAbsoluteError",
    "RootMeanSquaredError",
    "FalseDiscoveryRate",
    "FalseNegativeRate",
    "DetectionRate",
]




def average_precision_binary(labels: np.ndarray, scores: np.ndarray) -> float:
    """Compute non-interpolated average precision for binary labels.

    This avoids a hard dependency on scikit-learn. If no true positives exist in
    the truth vector, AUPRC is undefined and reported as NA.
    """
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=float)
    positive_count = int(labels.sum())
    if positive_count == 0 or len(labels) == 0:
        return np.nan

    finite_scores = np.where(np.isfinite(scores), scores, -np.inf)
    order = np.argsort(-finite_scores, kind="mergesort")
    sorted_labels = labels[order]
    sorted_scores = finite_scores[order]

    distinct = np.where(np.diff(sorted_scores))[0]
    threshold_idxs = np.r_[distinct, len(sorted_scores) - 1]
    true_positives = np.cumsum(sorted_labels)[threshold_idxs]
    false_positives = 1 + threshold_idxs - true_positives

    precision = true_positives / (true_positives + false_positives)
    recall = true_positives / positive_count
    return float(np.sum(np.diff(np.r_[0.0, recall]) * precision))


def safe_divide(numerator: float, denominator: float) -> float:
    if denominator == 0 or np.isclose(denominator, 0.0):
        return np.nan
    return float(numerator) / float(denominator)


def f1_score(precision: float, recall: float) -> float:
    if np.isnan(precision) or np.isnan(recall):
        return np.nan
    return safe_divide(2 * precision * recall, precision + recall)


def is_constant(values: np.ndarray) -> bool:
    return len(values) < 2 or np.allclose(values, values[0])


def compute_metrics(
    true_vec: pd.Series | np.ndarray,
    pred_vec: pd.Series | np.ndarray,
    presence_threshold: float = 0.0,
) -> dict[str, float | int]:
    """Compute all supported metrics on aligned vectors."""
    true_arr = np.asarray(true_vec, dtype=float)
    pred_arr = np.asarray(pred_vec, dtype=float)
    true_present = true_arr > presence_threshold
    pred_present = pred_arr > presence_threshold

    tp = int(np.logical_and(true_present, pred_present).sum())
    fp = int(np.logical_and(~true_present, pred_present).sum())
    fn = int(np.logical_and(true_present, ~pred_present).sum())
    tn = int(np.logical_and(~true_present, ~pred_present).sum())

    precision = safe_divide(tp, tp + fp)
    recall = safe_divide(tp, tp + fn)
    f1 = f1_score(precision, recall)
    fdr = safe_divide(fp, tp + fp)
    fnr = safe_divide(fn, tp + fn)
    detection_rate = recall
    accuracy = safe_divide(tp + tn, tp + fp + fn + tn)
    specificity = safe_divide(tn, tn + fp)
    negative_predictive_value = safe_divide(tn, tn + fn)
    mcc_denominator = float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    matthews = np.nan if np.isclose(mcc_denominator, 0.0) else float(((tp * tn) - (fp * fn)) / np.sqrt(mcc_denominator))
    auprc = average_precision_binary(true_present, pred_arr)

    overlap = float(np.minimum(true_arr, pred_arr).sum())
    true_sum = float(true_arr.sum())
    pred_sum = float(pred_arr.sum())
    max_sum = float(np.maximum(true_arr, pred_arr).sum())
    abundance_precision = safe_divide(overlap, pred_sum)
    abundance_recall = safe_divide(overlap, true_sum)
    abundance_f1 = f1_score(abundance_precision, abundance_recall)

    diff = true_arr - pred_arr
    abs_diff = np.abs(diff)
    l1 = float(abs_diff.sum())
    l2 = float(np.sqrt(np.square(diff).sum()))
    bc = 0.0 if np.isclose(true_sum + pred_sum, 0.0) else l1 / (true_sum + pred_sum)

    if np.isclose(true_sum, 0.0) and np.isclose(pred_sum, 0.0):
        rjsd = 0.0
    else:
        eps = 1e-12
        p = true_arr / true_sum if true_sum > 0 else np.zeros_like(true_arr)
        q = pred_arr / pred_sum if pred_sum > 0 else np.zeros_like(pred_arr)
        p = (p + eps) / (p + eps).sum()
        q = (q + eps) / (q + eps).sum()
        rjsd = float(jensenshannon(p, q))

    if is_constant(true_arr) or is_constant(pred_arr):
        spearman = np.nan
        pearson = np.nan
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            spearman = spearmanr(true_arr, pred_arr).correlation
            pearson_result = pearsonr(true_arr, pred_arr)
            pearson = getattr(pearson_result, "statistic", pearson_result[0])
        spearman = float(spearman) if np.isfinite(spearman) else np.nan
        pearson = float(pearson) if np.isfinite(pearson) else np.nan

    true_norm = float(np.linalg.norm(true_arr))
    pred_norm = float(np.linalg.norm(pred_arr))
    cosine = safe_divide(float(np.dot(true_arr, pred_arr)), true_norm * pred_norm)

    return {
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "Precision": precision,
        "Recall": recall,
        "F1": f1,
        "Spearman": spearman,
        "AbundancePrecision": abundance_precision,
        "AbundanceRecall": abundance_recall,
        "AbundanceF1": abundance_f1,
        "L1": l1,
        "L2": l2,
        "BC": bc,
        "rJSD": rjsd,
        "Jaccard": safe_divide(tp, tp + fp + fn),
        "WeightedJaccard": safe_divide(overlap, max_sum),
        "Pearson": pearson,
        "Cosine": cosine,
        "AUPRC": auprc,
        "Accuracy": accuracy,
        "Specificity": specificity,
        "NegativePredictiveValue": negative_predictive_value,
        "MatthewsCorrelationCoefficient": matthews,
        "MeanAbsoluteError": float(abs_diff.mean()) if len(abs_diff) else np.nan,
        "RootMeanSquaredError": float(np.sqrt(np.square(diff).mean())) if len(diff) else np.nan,
        "FalseDiscoveryRate": fdr,
        "FalseNegativeRate": fnr,
        "DetectionRate": detection_rate,
    }
