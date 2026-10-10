# Greengenes2 V3-V4 reference preparation — CLOSED
Date: 2026-10-06

Status:
ORIGINAL HISTORICAL CODE RECOVERED AND VERIFIED.

Original script:
scripts/12_gg2_extract_v3v4.pbs

Historical dated copy:
logs/scripts_snapshot_20260822_GG2_COMPLETE/12_gg2_extract_v3v4.pbs

The current script and dated historical snapshot are exact
byte-for-byte matches.

Reference extraction command:
qiime feature-classifier extract-reads

Input sequence reference:
10_greengenes_test/reference/2022.10.backbone.full-length.fna.qza

Primers:
338F = ACTCCTACGGGAGGCAGCAG
806R = GGACTACHVGGGTWTCTAAT

Extraction parameters:
- min length = 350
- max length = 550
- read orientation = both
- n_jobs = 1

Output reference:
10_greengenes_test/reference/gg2_2022_10_338F_806R.qza

Extraction statistics:
10_greengenes_test/reference/gg2_2022_10_338F_806R_stats.qza

Reference artifact SHA-256:

2022.10.backbone.full-length.fna.qza
07314c939f5a14162d889be57dd2813ea02ba80961a035df435d4273bc4a6893

2022.10.backbone.tax.qza
6495d6da4097ff51391ae7ce8df09c0ee00807c8bb64e4ce8cae1c9d5e07f5f0

gg2_2022_10_338F_806R.qza
0654c98c8b20cdaacd463eaf3b99ba49b629d80392879a9ae2bf97373fd1b851

gg2_2022_10_338F_806R_stats.qza
fd89478af55645c7652e0d095240b1ec962cde1a10ba5869575da50e3aa8e204
