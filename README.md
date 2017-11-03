# Regional_groundwater_chemistry_trends_in_R

This simple R script reads in selected analytes from the State of California’s groundwater Ambient Monitoring and Assessment (GAMA) database and performs principal component analysis and spatial correlation analysis on the resulting data subset. The script requires the following R packages:
* reshape
* corrplot
* ncf
The script requires data to be read in in one of two forms. The first consists of the raw, tab-delimited text files that are downloaded from the GAMA website (http://geotracker.waterboards.ca.gov/gama/datadownload). These are filtered by the script for only selected analytes. The filtered data set is then saved to a single pivoted, comma-delimited text file. This second file can then be read by the script in subsequent analyses, saving a lot of time. Neither the raw data files for the 59 counties, nor any composite data files, are supplied here because of the large file sizes. Therefore, the user will need to develop these as required.
Separately, a text file needs to be present that tells the script which county data sets to include in the analysis. An example, for the full set from all 59 California counties, is included with the script.
I’ve used this script to look in to the distributions of various major ion associations and trends in groundwater in California’s San Joaquin Valley, specifically with respect to nitrate. This application is described in a bit more detail on my blog (link is pending).
I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.
THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

