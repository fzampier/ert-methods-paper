| Scenario | Method | Inspect every | Crossing rate | Median crossing |
| --- | --- | --- | --- | --- |
| Alternative: 35% vs 30% | e-RTb adaptive |   1 | 45.9% (0.70pp) | 1,343 |
| Alternative: 35% vs 30% | e-RTb adaptive |  50 | 41.7% (0.70pp) | 1,525 |
| Alternative: 35% vs 30% | e-RTb adaptive | 100 | 40.0% (0.69pp) | 1,600 |
| Alternative: 35% vs 30% | e-RTb adaptive | 250 | 37.2% (0.68pp) | 1,750 |
| Alternative: 35% vs 30% | e-RTb design directional |   1 | 74.9% (0.61pp) | 1,282 |
| Alternative: 35% vs 30% | e-RTb design directional |  50 | 72.4% (0.63pp) | 1,350 |
| Alternative: 35% vs 30% | e-RTb design directional | 100 | 71.1% (0.64pp) | 1,400 |
| Alternative: 35% vs 30% | e-RTb design directional | 250 | 69.0% (0.65pp) | 1,500 |
| Alternative: 35% vs 30% | e-RTb mixture two-sided |   1 | 64.4% (0.68pp) | 1,506 |
| Alternative: 35% vs 30% | e-RTb mixture two-sided |  50 | 62.1% (0.69pp) | 1,600 |
| Alternative: 35% vs 30% | e-RTb mixture two-sided | 100 | 60.9% (0.69pp) | 1,600 |
| Alternative: 35% vs 30% | e-RTb mixture two-sided | 250 | 59.0% (0.70pp) | 1,750 |
| Null: 35% vs 35% | e-RTb adaptive |   1 | 3.1% (0.25pp) |   421 |
| Null: 35% vs 35% | e-RTb adaptive |  50 | 2.0% (0.20pp) |   625 |
| Null: 35% vs 35% | e-RTb adaptive | 100 | 1.5% (0.17pp) |   700 |
| Null: 35% vs 35% | e-RTb adaptive | 250 | 1.2% (0.15pp) |   750 |
| Null: 35% vs 35% | e-RTb design directional |   1 | 4.1% (0.28pp) | 1,358 |
| Null: 35% vs 35% | e-RTb design directional |  50 | 3.2% (0.25pp) | 1,400 |
| Null: 35% vs 35% | e-RTb design directional | 100 | 3.0% (0.24pp) | 1,500 |
| Null: 35% vs 35% | e-RTb design directional | 250 | 2.5% (0.22pp) | 1,500 |
| Null: 35% vs 35% | e-RTb mixture two-sided |   1 | 3.3% (0.25pp) | 1,655 |
| Null: 35% vs 35% | e-RTb mixture two-sided |  50 | 2.7% (0.23pp) | 1,750 |
| Null: 35% vs 35% | e-RTb mixture two-sided | 100 | 2.3% (0.21pp) | 1,750 |
| Null: 35% vs 35% | e-RTb mixture two-sided | 250 | 1.8% (0.19pp) | 1,750 |

Same simulated streams as the fair-comparator run (seeds 20260542/43); wealth computed patient-by-patient, then inspected only at multiples of the batch size (final look always included) -- batched INSPECTION of a sequentially updated process, at the point-in-time wealth. Inspecting a subset of times preserves anytime validity and is conservative under the null. Detection delay is NOT bounded by the batch size under this rule: an excursion that crosses 20 and falls back between inspections is missed permanently, not delayed (design wager: 74.9% patient-level vs 71.1% at every-100 inspection). If per-update wealth is logged and the inspector instead checks the running maximum since the last inspection, every patient-level crossing is caught at the next boundary (delay then bounded by the batch size, crossing set identical, validity unchanged -- Ville bounds the supremum). [Caption corrected 2026-08-13, audit; simulated numbers untouched.]
