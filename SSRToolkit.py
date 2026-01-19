#!/usr/bin/env python3
"""
SSRToolkit.py

- Detecta SSR perfectos en FASTA (solo dentro de bloques A/C/G/T; no cruza Ns)
- Extrae flancos limpios (se corta al primer caracter no-ACGT)
- Diseña primers flanqueando el SSR con primer3-py
- Reporte de:
    * SSR totales
    * SSR por tamaño de motivo (mono/di/tri/...)
    * SSR por motivo exacto (top N)
    * SSR evaluados para primers
    * SSR con primers diseñados
    * duplicados: primers left duplicados, right duplicados, pares (L+R) duplicados

Outputs:
  <prefix>.ssr_hits.tsv
  <prefix>.primers.tsv
  <prefix>.summary.txt

Requisitos:
  pip install biopython pandas primer3-py tqdm

Ejemplo:
  python SSRToolkit_no_multihit.py -i TAIR10.1_genomic.fna -o athaliana --threads 8 --flank 600
"""

import argparse
import re
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple, Any

import pandas as pd
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

try:
    import primer3
except ImportError:
    primer3 = None


ACGT_BLOCK_RE = re.compile(r"[ACGT]+", re.IGNORECASE)


@dataclass
class SSRHit:
    seq_id: str
    start_1based: int
    end_1based: int
    motif: str
    repeats: int
    ssr_len: int
    ssr_seq: str
    left_flank: str
    right_flank: str


def normalize_seq(seq: str) -> str:
    return seq.upper()


def parse_min_repeats_by_k(s: str) -> Dict[int, int]:
    """
    "1:10,2:6,3:5,4:5,5:4,6:4" -> {1:10,2:6,...}
    """
    out: Dict[int, int] = {}
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        k, v = part.split(":")
        out[int(k.strip())] = int(v.strip())
    return out


def get_clean_flanks(full_seq: str, start0: int, end0: int, flank: int) -> Tuple[str, str]:
    """
    Devuelve flancos A/C/G/T a cada lado sin cruzar Ns u otros símbolos.
    start0/end0 en coordenadas 0-based inclusive del SSR en full_seq.
    """
    n = len(full_seq)

    # left
    left_chars = []
    i = start0 - 1
    while i >= 0 and len(left_chars) < flank:
        c = full_seq[i]
        if c not in "ACGT":
            break
        left_chars.append(c)
        i -= 1
    left_flank = "".join(reversed(left_chars))

    # right
    right_chars = []
    i = end0 + 1
    while i < n and len(right_chars) < flank:
        c = full_seq[i]
        if c not in "ACGT":
            break
        right_chars.append(c)
        i += 1
    right_flank = "".join(right_chars)

    return left_flank, right_flank


def find_ssrs_perfect(
    full_seq: str,
    seq_id: str,
    flank: int,
    min_motif: int,
    max_motif: int,
    min_repeats_by_k: Dict[int, int],
    min_len: int,
    max_len: int,
    allow_overlap: bool = False,
) -> List[SSRHit]:
    """
    SSR perfectos dentro de bloques ACGT.
    """
    hits: List[SSRHit] = []
    full_seq = normalize_seq(full_seq)

    for block_match in ACGT_BLOCK_RE.finditer(full_seq):
        block_start0 = block_match.start()
        block_seq = block_match.group(0)
        block_len = len(block_seq)

        occupied = [False] * block_len  # evita solapes dentro del mismo bloque

        # Si min_repeats_by_k tiene keys definidas, solo buscar esos motivos
        motifs_to_search = list(min_repeats_by_k.keys()) if min_repeats_by_k else range(min_motif, max_motif + 1)
        
        for k in motifs_to_search:
            min_rep = min_repeats_by_k.get(k, 1)
            if min_rep < 1:
                continue

            pattern = re.compile(rf"(([ACGT]{{{k}}})\2{{{min_rep-1},}})")
            for m in pattern.finditer(block_seq):
                ssr_seq = m.group(1)
                motif = m.group(2)
                ssr_len = len(ssr_seq)

                if ssr_len < min_len or ssr_len > max_len:
                    continue

                start0_block = m.start(1)
                end0_block = m.end(1) - 1

                if not allow_overlap:
                    if any(occupied[i] for i in range(start0_block, end0_block + 1)):
                        continue
                    for i in range(start0_block, end0_block + 1):
                        occupied[i] = True

                start0_full = block_start0 + start0_block
                end0_full = block_start0 + end0_block

                repeats = ssr_len // k
                left_flank, right_flank = get_clean_flanks(full_seq, start0_full, end0_full, flank)

                hits.append(
                    SSRHit(
                        seq_id=seq_id,
                        start_1based=start0_full + 1,
                        end_1based=end0_full + 1,
                        motif=motif,
                        repeats=repeats,
                        ssr_len=ssr_len,
                        ssr_seq=ssr_seq,
                        left_flank=left_flank,
                        right_flank=right_flank,
                    )
                )

    hits.sort(key=lambda h: (h.seq_id, h.start_1based, -h.ssr_len))
    return hits


def design_primers_for_hit(
    hit: SSRHit,
    product_min: int,
    product_max: int,
    primer_min: int,
    primer_opt: int,
    primer_max: int,
    tm_min: float,
    tm_opt: float,
    tm_max: float,
    gc_min: float,
    gc_max: float,
    num_return: int = 1,
) -> Dict[str, Any]:
    """
    Diseña primers con primer3.
    Importante: incluye guardas para evitar crashes por template corto, etc.
    """
    if primer3 is None:
        raise RuntimeError("primer3-py no está instalado. Instala con: pip install primer3-py")

    template = hit.left_flank + hit.ssr_seq + hit.right_flank
    template_len = len(template)

    if template_len < product_min:
        return {
            "seq_id": hit.seq_id,
            "ssr_start": hit.start_1based,
            "ssr_end": hit.end_1based,
            "motif": hit.motif,
            "repeats": hit.repeats,
            "ssr_len": hit.ssr_len,
            "left_flank_len": len(hit.left_flank),
            "right_flank_len": len(hit.right_flank),
            "template_len": template_len,
            "primer_pair_found": 0,
            "skip_reason": f"template_len<{product_min}",
        }

    target_start = len(hit.left_flank)
    target_len = len(hit.ssr_seq)

    if target_len <= 0 or (target_start + target_len) > template_len:
        return {
            "seq_id": hit.seq_id,
            "ssr_start": hit.start_1based,
            "ssr_end": hit.end_1based,
            "motif": hit.motif,
            "repeats": hit.repeats,
            "ssr_len": hit.ssr_len,
            "left_flank_len": len(hit.left_flank),
            "right_flank_len": len(hit.right_flank),
            "template_len": template_len,
            "primer_pair_found": 0,
            "skip_reason": "invalid_target_region",
        }

    seq_args = {
        "SEQUENCE_ID": f"{hit.seq_id}:{hit.start_1based}-{hit.end_1based}",
        "SEQUENCE_TEMPLATE": template,
        "SEQUENCE_TARGET": [target_start, target_len],
    }

    global_args = {
        "PRIMER_OPT_SIZE": primer_opt,
        "PRIMER_MIN_SIZE": primer_min,
        "PRIMER_MAX_SIZE": primer_max,
        "PRIMER_OPT_TM": tm_opt,
        "PRIMER_MIN_TM": tm_min,
        "PRIMER_MAX_TM": tm_max,
        "PRIMER_MIN_GC": gc_min,
        "PRIMER_MAX_GC": gc_max,
        "PRIMER_PRODUCT_SIZE_RANGE": [[product_min, product_max]],
        "PRIMER_NUM_RETURN": int(num_return),
    }

    try:
        res = primer3.bindings.design_primers(seq_args, global_args)
    except OSError as e:
        return {
            "seq_id": hit.seq_id,
            "ssr_start": hit.start_1based,
            "ssr_end": hit.end_1based,
            "motif": hit.motif,
            "repeats": hit.repeats,
            "ssr_len": hit.ssr_len,
            "left_flank_len": len(hit.left_flank),
            "right_flank_len": len(hit.right_flank),
            "template_len": template_len,
            "primer_pair_found": 0,
            "skip_reason": f"primer3_oserror:{str(e)}",
        }
    except Exception as e:
        return {
            "seq_id": hit.seq_id,
            "ssr_start": hit.start_1based,
            "ssr_end": hit.end_1based,
            "motif": hit.motif,
            "repeats": hit.repeats,
            "ssr_len": hit.ssr_len,
            "left_flank_len": len(hit.left_flank),
            "right_flank_len": len(hit.right_flank),
            "template_len": template_len,
            "primer_pair_found": 0,
            "skip_reason": f"primer3_error:{type(e).__name__}:{str(e)}",
        }

    out = {
        "seq_id": hit.seq_id,
        "ssr_start": hit.start_1based,
        "ssr_end": hit.end_1based,
        "motif": hit.motif,
        "repeats": hit.repeats,
        "ssr_len": hit.ssr_len,
        "left_flank_len": len(hit.left_flank),
        "right_flank_len": len(hit.right_flank),
        "template_len": template_len,
        "primer_pair_found": int(res.get("PRIMER_PAIR_NUM_RETURNED", 0) > 0),
        "skip_reason": "",
    }

    if out["primer_pair_found"]:
        out.update(
            {
                "primer_left_seq": res["PRIMER_LEFT_0_SEQUENCE"],
                "primer_right_seq": res["PRIMER_RIGHT_0_SEQUENCE"],
                "primer_left_tm": res["PRIMER_LEFT_0_TM"],
                "primer_right_tm": res["PRIMER_RIGHT_0_TM"],
                "primer_left_gc": res["PRIMER_LEFT_0_GC_PERCENT"],
                "primer_right_gc": res["PRIMER_RIGHT_0_GC_PERCENT"],
                "product_size": res["PRIMER_PAIR_0_PRODUCT_SIZE"],
            }
        )

    return out


def process_record(job):
    """
    Worker para paralelizar por contig:
    - SSR mining
    - primer3 por SSR
    """
    seq_id, seq, params = job

    hits = find_ssrs_perfect(
        full_seq=seq,
        seq_id=seq_id,
        flank=params["flank"],
        min_motif=params["min_motif"],
        max_motif=params["max_motif"],
        min_repeats_by_k=params["min_repeats_by_k"],
        min_len=params["min_len"],
        max_len=params["max_len"],
        allow_overlap=params["allow_overlap"],
    )
    hit_dicts = [asdict(h) for h in hits]

    primer_rows: List[Dict[str, Any]] = []
    if not params["no_primers"] and hits:
        for h in hits:
            # flancos mínimos para que quepa un primer de largo primer_max
            if len(h.left_flank) < params["primer_max"] or len(h.right_flank) < params["primer_max"]:
                primer_rows.append(
                    {
                        "seq_id": h.seq_id,
                        "ssr_start": h.start_1based,
                        "ssr_end": h.end_1based,
                        "motif": h.motif,
                        "repeats": h.repeats,
                        "ssr_len": h.ssr_len,
                        "left_flank_len": len(h.left_flank),
                        "right_flank_len": len(h.right_flank),
                        "template_len": len(h.left_flank) + len(h.ssr_seq) + len(h.right_flank),
                        "primer_pair_found": 0,
                        "skip_reason": f"flank_len<{params['primer_max']}",
                    }
                )
                continue

            row = design_primers_for_hit(
                h,
                product_min=params["product_min"],
                product_max=params["product_max"],
                primer_min=params["primer_min"],
                primer_opt=params["primer_opt"],
                primer_max=params["primer_max"],
                tm_min=params["tm_min"],
                tm_opt=params["tm_opt"],
                tm_max=params["tm_max"],
                gc_min=params["gc_min"],
                gc_max=params["gc_max"],
                num_return=1,
            )
            primer_rows.append(row)

    return hit_dicts, primer_rows


def safe_tqdm(iterable, **kwargs):
    if tqdm is None:
        return iterable
    return tqdm(iterable, **kwargs)


def main():
    ap = argparse.ArgumentParser(description="SSR mining + Primer3 (sin minimap2) + reporte")
    ap.add_argument("-i", "--input", required=True, help="FASTA (genoma/contigs)")
    ap.add_argument("-o", "--out-prefix", default="ssr_out", help="Prefijo de salida")
    ap.add_argument("--threads", type=int, default=1, help="Procesos para paralelizar por contig (default 1)")

    ap.add_argument("--flank", type=int, default=200, help="bp max a cada lado del SSR (se corta si hay N)")
    ap.add_argument("--min-motif", type=int, default=1)
    ap.add_argument("--max-motif", type=int, default=6)
    ap.add_argument("--min-repeats-by-k", default="1:10,2:6,3:5,4:5,5:4,6:4")
    ap.add_argument("--min-len", type=int, default=10)
    ap.add_argument("--max-len", type=int, default=200)
    ap.add_argument("--allow-overlap", action="store_true")

    ap.add_argument("--no-primers", action="store_true")
    ap.add_argument("--product-min", type=int, default=100)
    ap.add_argument("--product-max", type=int, default=300)
    ap.add_argument("--primer-min", type=int, default=18)
    ap.add_argument("--primer-opt", type=int, default=20)
    ap.add_argument("--primer-max", type=int, default=24)
    ap.add_argument("--tm-min", type=float, default=57.0)
    ap.add_argument("--tm-opt", type=float, default=60.0)
    ap.add_argument("--tm-max", type=float, default=63.0)
    ap.add_argument("--gc-min", type=float, default=40.0)
    ap.add_argument("--gc-max", type=float, default=60.0)

    ap.add_argument("--top-motifs", type=int, default=30, help="Top N motivos a reportar (default 30)")

    args = ap.parse_args()

    if not args.no_primers and primer3 is None:
        raise RuntimeError("Falta primer3-py. Instala con: pip install primer3-py")

    min_repeats_by_k = parse_min_repeats_by_k(args.min_repeats_by_k)
    
    # Si se especifica --min-repeats-by-k con valores concretos, ajustar min_motif y max_motif
    min_motif = args.min_motif
    max_motif = args.max_motif
    if min_repeats_by_k:  # Si el usuario especificó valores concretos
        motif_keys = list(min_repeats_by_k.keys())
        min_motif = min(motif_keys)
        max_motif = max(motif_keys)

    params = {
        "flank": args.flank,
        "min_motif": min_motif,
        "max_motif": max_motif,
        "min_repeats_by_k": min_repeats_by_k,
        "min_len": args.min_len,
        "max_len": args.max_len,
        "allow_overlap": args.allow_overlap,
        "no_primers": args.no_primers,
        "product_min": args.product_min,
        "product_max": args.product_max,
        "primer_min": args.primer_min,
        "primer_opt": args.primer_opt,
        "primer_max": args.primer_max,
        "tm_min": args.tm_min,
        "tm_opt": args.tm_opt,
        "tm_max": args.tm_max,
        "gc_min": args.gc_min,
        "gc_max": args.gc_max,
    }

    # Preparar jobs por contig
    jobs = []
    for rec in SeqIO.parse(args.input, "fasta"):
        seq = normalize_seq(str(rec.seq))
        jobs.append((rec.id, seq, params))

    # Ejecutar (SSR + primers)
    all_hit_dicts: List[Dict[str, Any]] = []
    all_primer_rows: List[Dict[str, Any]] = []

    if args.threads <= 1:
        for j in safe_tqdm(jobs, desc="Procesando contigs (SSR+Primer3)", unit="contig"):
            hd, pr = process_record(j)
            all_hit_dicts.extend(hd)
            all_primer_rows.extend(pr)
    else:
        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            futures = [ex.submit(process_record, j) for j in jobs]
            for fu in safe_tqdm(as_completed(futures), total=len(futures), desc="Procesando contigs (SSR+Primer3)", unit="contig"):
                hd, pr = fu.result()
                all_hit_dicts.extend(hd)
                all_primer_rows.extend(pr)

    # Guardar outputs
    ssr_df = pd.DataFrame(all_hit_dicts)
    primers_df = pd.DataFrame(all_primer_rows)

    # Generar IDs únicos para primers
    if not primers_df.empty:
        primers_df.insert(0, "primer_id", [f"PRIMER_{i:06d}" for i in range(1, len(primers_df) + 1)])

    ssr_path = f"{args.out_prefix}.ssr_hits.tsv"
    primers_path = f"{args.out_prefix}.primers.tsv"
    duplicates_path = f"{args.out_prefix}.duplicates.tsv"
    summary_path = f"{args.out_prefix}.summary.txt"

    ssr_df.to_csv(ssr_path, sep="\t", index=False)
    primers_df.to_csv(primers_path, sep="\t", index=False)

    # ---------------------------
    # Reporte / métricas
    # ---------------------------
    total_ssr = int(len(ssr_df))

    # Por tamaño del motivo (1..6 etc.)
    by_k = {}
    by_motif = {}
    if total_ssr > 0:
        ssr_df["motif_len"] = ssr_df["motif"].astype(str).str.len()
        by_k = ssr_df["motif_len"].value_counts().sort_index().to_dict()
        by_motif = ssr_df["motif"].value_counts().to_dict()

    total_rows_primers = int(len(primers_df))
    primers_attempted = int(total_rows_primers)  # cada SSR genera una fila en primers_df en este diseño
    primers_ok = int((primers_df.get("primer_pair_found", 0) == 1).sum()) if total_rows_primers else 0

    # Duplicados (solo entre los OK)
    left_dup_total = right_dup_total = pair_dup_total = 0
    left_dup_unique_seqs = right_dup_unique_seqs = pair_dup_unique_seqs = 0

    duplicates_list = []

    if total_rows_primers and "primer_left_seq" in primers_df.columns:
        ok_df = primers_df[primers_df["primer_pair_found"] == 1].copy()

        if not ok_df.empty:
            # Left duplicates
            left_counts = ok_df["primer_left_seq"].value_counts()
            left_dup_unique_seqs = int((left_counts > 1).sum())
            left_dup_total = int(left_counts[left_counts > 1].sum())

            # Agregar left duplicates a la lista
            for left_seq in left_counts[left_counts > 1].index:
                dup_rows = ok_df[ok_df["primer_left_seq"] == left_seq]
                for _, row in dup_rows.iterrows():
                    duplicates_list.append({
                        "primer_id": row.get("primer_id", ""),
                        "seq_id": row["seq_id"],
                        "ssr_start": row["ssr_start"],
                        "ssr_end": row["ssr_end"],
                        "duplicate_type": "LEFT",
                        "primer_left_seq": left_seq,
                        "primer_right_seq": row.get("primer_right_seq", ""),
                        "motif": row["motif"],
                        "repeats": row["repeats"],
                    })

            # Right duplicates
            right_counts = ok_df["primer_right_seq"].value_counts()
            right_dup_unique_seqs = int((right_counts > 1).sum())
            right_dup_total = int(right_counts[right_counts > 1].sum())

            # Agregar right duplicates a la lista
            for right_seq in right_counts[right_counts > 1].index:
                dup_rows = ok_df[ok_df["primer_right_seq"] == right_seq]
                for _, row in dup_rows.iterrows():
                    duplicates_list.append({
                        "primer_id": row.get("primer_id", ""),
                        "seq_id": row["seq_id"],
                        "ssr_start": row["ssr_start"],
                        "ssr_end": row["ssr_end"],
                        "duplicate_type": "RIGHT",
                        "primer_left_seq": row.get("primer_left_seq", ""),
                        "primer_right_seq": right_seq,
                        "motif": row["motif"],
                        "repeats": row["repeats"],
                    })

            # Pair duplicates (L+R exactos)
            ok_df["primer_pair"] = ok_df["primer_left_seq"].astype(str) + "||" + ok_df["primer_right_seq"].astype(str)
            pair_counts = ok_df["primer_pair"].value_counts()
            pair_dup_unique_seqs = int((pair_counts > 1).sum())
            pair_dup_total = int(pair_counts[pair_counts > 1].sum())

            # Agregar pair duplicates a la lista
            for pair_seq in pair_counts[pair_counts > 1].index:
                left_seq, right_seq = pair_seq.split("||")
                dup_rows = ok_df[ok_df["primer_pair"] == pair_seq]
                for _, row in dup_rows.iterrows():
                    duplicates_list.append({
                        "primer_id": row.get("primer_id", ""),
                        "seq_id": row["seq_id"],
                        "ssr_start": row["ssr_start"],
                        "ssr_end": row["ssr_end"],
                        "duplicate_type": "PAIR",
                        "primer_left_seq": left_seq,
                        "primer_right_seq": right_seq,
                        "motif": row["motif"],
                        "repeats": row["repeats"],
                    })

    # Guardar duplicados en archivo
    if duplicates_list:
        duplicates_df = pd.DataFrame(duplicates_list)
        # Eliminar duplicados exactos (misma fila puede aparecer en múltiples categorías)
        duplicates_df = duplicates_df.drop_duplicates(subset=["primer_id", "duplicate_type"])
        duplicates_df.to_csv(duplicates_path, sep="\t", index=False)
    else:
        # Crear archivo vacío o con solo headers
        empty_dup_df = pd.DataFrame(columns=["primer_id", "seq_id", "ssr_start", "ssr_end", "duplicate_type", "primer_left_seq", "primer_right_seq", "motif", "repeats"])
        empty_dup_df.to_csv(duplicates_path, sep="\t", index=False)

    # Generar archivos de primers únicos y duplicados
    unique_primers_path = f"{args.out_prefix}.primers_unique.tsv"
    duplicated_primers_path = f"{args.out_prefix}.primers_duplicated.tsv"
    
    if total_rows_primers and "primer_left_seq" in primers_df.columns:
        # Obtener IDs de primers con duplicados
        dup_primer_ids = set(duplicates_df['primer_id'].unique()) if duplicates_list else set()
        
        # Filtrar solo primers exitosos
        successful_primers = primers_df[primers_df['primer_pair_found'] == 1].copy()
        
        # Separar únicos y duplicados
        unique_primers = successful_primers[~successful_primers['primer_id'].isin(dup_primer_ids)]
        duplicated_primers = successful_primers[successful_primers['primer_id'].isin(dup_primer_ids)]
        
        # Guardar archivos
        unique_primers.to_csv(unique_primers_path, sep="\t", index=False)
        duplicated_primers.to_csv(duplicated_primers_path, sep="\t", index=False)
    else:
        # Crear archivos vacíos si no hay primers
        empty_df = pd.DataFrame()
        empty_df.to_csv(unique_primers_path, sep="\t", index=False)
        empty_df.to_csv(duplicated_primers_path, sep="\t", index=False)

    # Escribir summary.txt
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("=== SSRToolkit Reporte (sin minimap2) ===\n\n")

        f.write(f"Input FASTA: {args.input}\n")
        f.write(f"Output prefix: {args.out_prefix}\n\n")

        f.write("---- SSR ----\n")
        f.write(f"SSR totales encontrados: {total_ssr}\n")
        f.write("SSR por tamaño de motivo (k):\n")
        if by_k:
            for k in sorted(by_k.keys()):
                name = {1:"mono",2:"di",3:"tri",4:"tetra",5:"penta",6:"hexa"}.get(k, f"k={k}")
                f.write(f"  - {name} ({k}): {by_k[k]}\n")
        else:
            f.write("  (sin SSR)\n")

        f.write(f"\nTop {args.top_motifs} motivos (motif exacto):\n")
        if by_motif:
            top = sorted(by_motif.items(), key=lambda x: x[1], reverse=True)[: args.top_motifs]
            for motif, c in top:
                f.write(f"  - {motif}: {c}\n")
        else:
            f.write("  (sin SSR)\n")

        f.write("\n---- Primers (Primer3) ----\n")
        f.write(f"Filas en primers.tsv: {total_rows_primers}\n")
        f.write(f"SSR evaluados para primers: {primers_attempted}\n")
        f.write(f"SSR con primers diseñados (primer_pair_found==1): {primers_ok}\n")
        if primers_attempted > 0:
            f.write(f"Éxito (primers_ok / evaluados): {100*primers_ok/primers_attempted:.2f}%\n")
        if total_ssr > 0:
            f.write(f"Éxito global (primers_ok / SSR totales): {100*primers_ok/total_ssr:.2f}%\n")

        f.write("\n---- Duplicados (solo entre primers OK) ----\n")
        f.write(f"Left primer: secuencias duplicadas (distintas): {left_dup_unique_seqs}\n")
        f.write(f"Left primer: ocurrencias totales dentro de duplicados: {left_dup_total}\n")
        f.write(f"Right primer: secuencias duplicadas (distintas): {right_dup_unique_seqs}\n")
        f.write(f"Right primer: ocurrencias totales dentro de duplicados: {right_dup_total}\n")
        f.write(f"Par (L+R): pares duplicados (distintos): {pair_dup_unique_seqs}\n")
        f.write(f"Par (L+R): ocurrencias totales dentro de duplicados: {pair_dup_total}\n")

        f.write("\nArchivos:\n")
        f.write(f"  - {ssr_path}\n")
        f.write(f"  - {primers_path}\n")
        f.write(f"  - {unique_primers_path} (Primers únicos sin duplicados)\n")
        f.write(f"  - {duplicated_primers_path} (Primers con duplicados)\n")
        f.write(f"  - {duplicates_path} (Listado detallado de duplicados)\n")
        f.write(f"  - {summary_path}\n")

    # Print corto
    print("\n=== RESUMEN ===")
    print(f"SSR totales: {total_ssr}")
    print(f"Primers OK: {primers_ok} / evaluados {primers_attempted}")
    print(f"Duplicados pares (L+R) distintos: {pair_dup_unique_seqs} (ocurrencias: {pair_dup_total})")
    
    # Estadísticas de únicos vs duplicados
    if total_rows_primers:
        successful = (primers_df.get("primer_pair_found", 0) == 1).sum()
        if "primer_id" in primers_df.columns:
            dup_ids = len(set(duplicates_df['primer_id'].unique())) if duplicates_list else 0
            unique_count = successful - dup_ids
            if unique_count > 0:
                print(f"\nPrimers únicos: {unique_count} ({100*unique_count/successful:.2f}%)")
            if dup_ids > 0:
                print(f"Primers duplicados: {dup_ids} ({100*dup_ids/successful:.2f}%)")
    
    print(f"\n✓ Archivos generados:")
    print(f"  - {ssr_path}")
    print(f"  - {primers_path}")
    print(f"  - {unique_primers_path}")
    print(f"  - {duplicated_primers_path}")
    print(f"  - {duplicates_path}")
    print(f"\nReporte: {summary_path}")


if __name__ == "__main__":
    main()
