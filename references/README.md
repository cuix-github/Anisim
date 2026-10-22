# PF-FLIP and ST-FLIP literature library

The core papers and the six references used in PF-FLIP Section 7 are grouped here by their role in the implementation. The [October 22 audit](../docs/literature/spray-adaptivity-audit.md) explains the distinction between spatial adaptivity, spray dynamics and render-time reconstruction, and what the existing 40-second cache can support.

**As of 2026-10-22: five complete paper PDFs, one authorized book preview, and two entries without a local PDF.** Downloaded does not mean fully read or fully implemented. The [catalog](catalog.json) records exact versions, authors, source URLs, SHA-256 hashes, PDF page counts and the actual reading scope. No source is marked as exhaustively read in this audit.

The public repository contains this index, metadata and analysis. Reading copies live under the ignored `references/papers/` directory; a fresh clone will therefore need the linked author/publisher downloads. Public availability of a PDF is not treated as permission to republish it. Existing simulation caches, images and scripts stay where they were.

```text
references/
├── README.md
├── catalog.json
└── papers/                              local, excluded from Git
    ├── 01_core/
    │   ├── Adaptive_Phase_Field_FLIP_preprint.pdf
    │   └── ST-FLIP.pdf
    ├── 02_spray_and_droplets/
    │   ├── Gissler_2017_Generalized_Drag_Force.pdf
    │   ├── Ihmsen_2012_Unified_Spray_Foam_Bubbles.pdf
    │   ├── Sirignano_2010_Droplets_Sprays_PREVIEW_ONLY.pdf
    │   ├── song_2005_SOURCE.txt
    │   ├── babinsky_sojka_2002_SOURCE.txt
    │   ├── sirignano_2010_SOURCE.txt
    │   └── *_access_status.json
    └── 03_surface_reconstruction/
        └── Bhattacharya_2011_Level_Set_Skinning.pdf
```

## Core methods

| Paper | Local reading copy | What to read for this question |
| --- | --- | --- |
| Braun, Bender, Thuerey (2025), [Adaptive Phase-Field-FLIP for Very Large Scale Two-Phase Fluid Simulation](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf) | [Full author PDF, 23 pages](papers/01_core/Adaptive_Phase_Field_FLIP_preprint.pdf) | §3.4 escape; §5–6 spatial adaptivity; §7 spray and rendering; Table 7. Targeted sections reread, not a complete reread |
| Braun, Winchenbach, Bender, Thuerey (2026), [Spatiotemporal FLIP for Fast Free-Surface and Two-Phase Simulation With Very Large Time Steps](https://ge.in.tum.de/download/ST-FLIP.pdf) | [Full author PDF, 20 pages](papers/01_core/ST-FLIP.pdf) | p.9 implementation/time synchronization; Appendix B water reconstruction. Targeted reading |

## Spray and droplets

| Reference | Local status | Role and actual reading scope |
| --- | --- | --- |
| Gissler et al. (2017), [Generalized drag force for particle-based simulations](https://cg.informatik.uni-freiburg.de/publications/2017_CAG_generalizedDragForce_v2.pdf) | [Full author PDF, 12 pages](papers/02_spray_and_droplets/Gissler_2017_Generalized_Drag_Force.pdf) | Relative air velocity and drag. Abstract and §3–4.2, especially author PDF pp.3–5; not the entire paper |
| Song, Shin, Ko (2005), [Stable but nondissipative water](https://snu.elsevierpure.com/en/publications/stable-but-nondissipative-water/) | No local PDF; [source pointer](papers/02_spray_and_droplets/song_2005_SOURCE.txt) | Lagrangian fragment/droplet precedent. Abstract and §5.1–5.4 in [author-uploaded web text](https://www.researchgate.net/publication/220183677_Stable_but_nondissipative_water); PDF endpoint returned HTTP 403 |
| Ihmsen et al. (2012), [Unified spray, foam and air bubbles for particle-based fluids](https://cg.informatik.uni-freiburg.de/publications/2012_CGI_sprayFoamBubbles.pdf) | [Full author PDF, 9 pages](papers/02_spray_and_droplets/Ihmsen_2012_Unified_Spray_Foam_Bubbles.pdf) | Secondary-particle generation and single-direction coupling. Abstract, introduction and §3.1–3.2. Author PDF title omits “air”; publisher title includes it |
| Babinsky and Sojka (2002), [Modeling drop size distributions](https://www.sciencedirect.com/science/article/abs/pii/S0360128502000047) | No local PDF; [source pointer](papers/02_spray_and_droplets/babinsky_sojka_2002_SOURCE.txt) | Radius-distribution background. Publisher abstract and introduction only; no legitimate public full PDF found |
| Sirignano (2010), [Fluid Dynamics and Transport of Droplets and Sprays, second edition](https://www.cambridge.org/core/books/fluid-dynamics-and-transport-of-droplets-and-sprays/EEEF907721AE679B11BC62AA59C61DE2) | [Authorized preview only, 49 PDF pages](papers/02_spray_and_droplets/Sirignano_2010_Droplets_Sprays_PREVIEW_ONLY.pdf) | Evaporation. Targeted printed pp.21–24, especially Eq.2.21. Preview includes frontmatter and printed pp.1–29; the full book was not downloaded or read |

PF-FLIP cites components of these methods. That does not establish adoption of every force, emitter, feedback rule or classification used by the referenced authors. In particular, Ihmsen's offline secondary particles demonstrate a possible workflow; they do not prove that any PF-FLIP cache contains enough information for exact air-driven spray replay.

## Water surface reconstruction

| Reference | Local status | Role and actual reading scope |
| --- | --- | --- |
| Bhatacharya, Gao, Bargteil (2011), [A Level-set Method for Skinning Animated Particle Data](https://cal.cs.umbc.edu/Papers/Bhattacharya-2011-ALM/Bhattacharya-2011-ALM.pdf) | [Full SCA 2011 author PDF, 8 pages](papers/03_surface_reconstruction/Bhattacharya_2011_Level_Set_Skinning.pdf) | Bounded level-set skinning and smoothing. Abstract, §3.1 and opening §3.2. The title page spells “Bhatacharya”; the institutional page uses “Bhattacharya.” This is the cited 2011 paper, not its 2015 expansion |

## Reading order and correspondence checks

Start with PF-FLIP §7 alongside its escape definition in §3.4, then inspect Ihmsen's birth sampling, Gissler's drag inputs and Sirignano's evaporation assumptions. Read PF-FLIP §5–6 separately to evaluate the memory/performance architecture. Use ST-FLIP Appendix B and the 2011 skinning paper for the water surface. This keeps optical representation, particle dynamics and pressure-grid refinement from being mistaken for one operation.

For each implemented component, record the paper's requirement, the code path, the actual production configuration, a validation result and the cache fields required to reproduce it. The current [audit](../docs/literature/spray-adaptivity-audit.md) begins this mapping. A missing detail in the paper should stay an author question or a labeled local design choice, not an undocumented claim of exact reproduction.
