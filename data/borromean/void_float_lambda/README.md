# Void: step 2's first scan, run from commit e94635f

These 11 files are **not results**.

**The cause.** `step2_matched.py prep` stored Λ(p^k) as float64 numbers. The prime terms therefore carried absolute errors of about 1e-17, and the eigenvalues, near e^{−150} or smaller, sat on a noise floor. The values are around 1e-18 with random signs.

**Why the precision loop missed it.** `stable_lams` raises the working precision but cannot see input errors that do not change with precision, so it accepted the noise as stable.

**How it was caught.** Negative minima appeared at several points for both objects.

**The fix.** It is in commit 8422f34:
- integer coefficients c with Λ(p^k) = c·log p, with log p formed at working precision;
- gates V1 and V2 on the production path, with a float64 positive control that must fire (`data/borromean/step2_validate.json`).

**The registration.** The predictions in docs/BORROMEAN.md were not changed. These values carry no information about λ.
