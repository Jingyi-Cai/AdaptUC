# AdaptUC

**AdaptUC** is a computational framework for designing microbial strains capable of assimilating previously unadapted carbon sources (e.g., methanol, formate, xylose) through Adaptive Laboratory Evolution (ALE). It leverages genome-scale metabolic models and bi-level mixed-integer programming (MIP) to predict gene knockout strategies that maximize the evolutionary driving force (SADF) for substrate adaptation.

---

## Dependencies

- **MATLAB** (tested with R2019b)
- **COBRA Toolbox** (tested with version 2.42.0)
- **IBM CPLEX** solver (tested with version 12.10)

Ensure all dependencies are installed and accessible in your MATLAB environment before running AdaptUC.

---

##  Quick Start

1. **Clone the repository**:
   ```bash
   git clone https://github.com/Jingyi-Cai/AdaptUC.git
   cd AdaptUC
   ```

2. **Open MATLAB** and navigate to the project folder.

3. **Run the case study script**:
   ```matlab
   AdaptUC_case_study.m
   ```
   This script will:
   - Load the *E. coli* iML1515 model (`iML1515.mat`).
   - Integrate methanol assimilation pathway reactions.
   - Define environmental conditions and knockout candidate set (`delCand.mat`).
   - Configure AdaptUC options (maximum knockouts, growth ratios, etc.).
   - Solve the bi-level MIP and print the optimal knockout strategies.

---
## License

AdaptUC is licensed under the Apache License, Version 2.0. You may obtain a copy at:

```
http://www.apache.org/licenses/LICENSE-2.0
```

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
## Input Files

- `iML1515.mat` : Genome-scale model of *E. coli* K-12 (iML1515)
- `delCand.mat` : Pre-filtered list of reaction deletion candidates
- `AdaptUC_case_study.m` : Main MATLAB script for running the case study

---


## Parameter Configuration

In `AdaptUC_case_study.m`, adjust `optionstcf` fields to explore different scenarios:

- `optionstcf.maxdel` : Maximum number of knockouts 
- `optionstcf.c1lb` : Minimum growth ratio on the unadapted substrate
- `optionstcf.co_sub_growthRatio` : Maximum reduced growth ratio on co-substrate
- `optionstcf.co_sub_improve` : Required fold-improvement with co-substrate present
- *(See script comments for full list)*
