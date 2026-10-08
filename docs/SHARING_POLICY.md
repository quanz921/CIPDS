# Sharing policy

No actual hospital or NHANES participant records are distributed. The generator reads no patient records, samples none, and is not fitted to real data or empirical distributions. It uses hand-written constants and a seed.

Excluded: raw data, linkage keys, dates, individual predictions, OOF tables, frozen hospital models, serialized training objects, logs, manuscripts and Git history. Private GitHub status does not replace these safeguards.

Synthetic outputs demonstrate selected computations, not the paper's findings or independent numerical reproduction. Exact reproduction needs authorised institutional inputs and model artifacts in a controlled environment. NHANES is public but is not bundled. NHANES alone cannot reproduce hospital-trained scores without the restricted fitted models or their authorised regeneration.

TRIPOD+AI separates data and code sharing. NIH guidance recognises disclosure risk after de-identification and possible controlled access. These sources inform the approach; no NIH applicability or privacy certification is claimed.

- [TRIPOD+AI](https://www.bmj.com/content/385/bmj-2023-078378)
- [NIH privacy guidance](https://www.grants.nih.gov/grants/guide/notice-files/NOT-OD-22-213.html)
- [NHANES](https://www.cdc.gov/nchs/nhanes/)
- [Linked mortality](https://www.cdc.gov/nchs/data-linkage/mortality-public.htm)

Checks include an allowlist, synthetic markers, hashes, syntax checks and credential/path scans. They are packaging checks, not a differential-privacy proof.
