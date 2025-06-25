# Plants in limbo —— only macroalgal detritus can remain viable for months
Repository accompanying Research in Context submission to special issue Plant Senescence in *Annals of Botany*

This reporsitory contains two folders, one for the research component (`Seagrass`) and one for the context component (`Meta-analysis`):

1. `Seagrass`
   - `Oxygen` Folder of csv files exported from the oxygen meter
   - `Prior`
      - `Screenshots` Folder of plot-digitised figures
      - `Prior.csv`
        - `Reference` Citation
        - `DOI` Digital object identifier
        - `Species` "*Halophila ovalis*" or "*Amphibolis antarctica*"
        - `Variable` "Light-saturated net photosynthesis" or "Detrital respiration"
        - `Flux` Measurement given in µmol O<sub>2</sub> g<sup>-1</sup> fresh mass h<sup>-1</sup>
    - `Mass.csv`
      - `Date` Date given as dd.mm.yy
      - `Species` "*Halophila ovalis*" or "*Amphibolis antarctica*"
      - `Leaf` Leaf replicate
      - `Mass` Blotted mass given in g
    - `PAR_Incubation.csv`
      - `PAR` Photosynthetically active radiation during O<sub>2</sub> measurement given in µmol photons m<sup>-2</sup> s<sup>-1</sup>
    - `PAR_Tank.csv`
      - `PAR` Photosynthetically active radiation in experimental tank given in µmol photons m<sup>-2</sup> s<sup>-1</sup>
    - `Volume.csv`
      - `Volume` Incubation jar volume given in L, determined by weighing ultrapure water
    - `Seagrass.R` R script to analyse all seagrass data
3. `Meta-analysis`
   - `Screenshots` Folder of plot-digitised figures
   - `Meta.csv` Detrital photosynthesis meta-dataset (see https://github.com/lukaseamus/detrital-photosynthesis)
   - `Meta.R` R script to perform the meta-analysis

Luka Seamus Wright, 20th March 2025
