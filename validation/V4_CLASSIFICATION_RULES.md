# V4 classification contract

This document implements V0 using archived evidence. Keys are versioned Ensembl gene ID, model and coefficient/test. Panel membership does not duplicate keys in the main table.

Precedence is: (1) missing required evidence or unavailable/nonsignificant baseline -> insufficient information / not estimable; (2) baseline-significant with any applicable recorded influence/addback/weighting effect/sign flag or FDR crossing -> sensitivity-qualified; (3) otherwise -> supported across the specified checks. Nonsignificant is not synonymous with computational non-estimability and is not proof of a null effect.

Primary models require V2 plus seven V3 cases. Conditional/restricted models use V2 and matching-cohort/formula weighting evidence; V3 was not in their frozen scope and is explicitly N/A. Subtype/combined adjustment changes the estimand and contextualizes rather than automatically vetoing a primary association. FDR boundary sensitivity alone is not proof of biological instability. No missing/out-of-scope test is treated as a pass. The three non-I primary Stage contrasts now have both V3 and verified V2 supplement evidence. Their classification follows the same precedence; baseline nonsignificance and recorded sensitivity flags remain consequential.

Stage omnibus is a joint F test without coefficient direction. Individual Stage contrasts use their own saved moderated t-test FDR. No omnibus FDR is copied to individual contrasts. Every FDR value is read from an archived result; V4 performs no new hypothesis test or BH computation.

Amendment note: decision rules and thresholds are unchanged. Only the historical coverage statement above is updated. The original contract is retained as V4_CLASSIFICATION_RULES_original.md.
