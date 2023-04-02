#include<fstream>
#include<iostream>
#include<string>
#include<vector>
#include<cstdlib>
#include<math.h>
#define PI 3.1415926535897932384626433832795
#define ZERO 1e-6
using namespace std;

vector<string> parseLine(string input, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = input.find(delimiter, pos_start)) != std::string::npos) {
        token = input.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (input.substr (pos_start));
    return res;
}
struct point {
    int x;
    int y;
    point(int x, int y) {
        this->x = x;
        this->y = y;
    }
};

class Args {
    public:
        int Imax = 128;
        int Jmax = 128;
        float rho = 1.0;
        float lineSpacing = 1.0;
        float eps = 1e-10;
        int numOfAngles = 1000;
        int numOfIterations = 10;
        float dx = 0.006;
        float dy = 0.006;
        string dataFile = "dataFile.txt";
        int numOfHeaderLines = 3;
        string outFile = "out.txt";
        bool hasRealP = false;       

    public:
        Args() {
            
        }

        void setFieldFromPair(vector<string> pair) {
            string key = pair[0];
            if (key == "Imax") {
                Imax = std::stoi(pair[1]);
            } else if (key == "Jmax") {
                Jmax = std::stoi(pair[1]);
            } else if (key == "rho") {
                rho =  std::stof(pair[1]);
            } else if (key == "lineSpacing") {
                lineSpacing = std::stof(pair[1]);
            } else if (key == "eps") {
                eps = std::stof(pair[1]);
            } else if (key == "numOfAngles") {
                numOfAngles = std::stoi(pair[1]);
            } else if (key == "numOfIterations") {
                numOfIterations = std::stoi(pair[1]);
            } else if (key == "dx") {
                dx = std::stof(pair[1]);
            } else if (key == "dy") {
                dy = std::stof(pair[1]);
            } else if (key == "dataFile") {
                dataFile = pair[1];
            } else if (key == "numOfHeaderLines") {
                numOfHeaderLines = std::stoi(pair[1]);
            } else if (key == "outFile") {
                outFile = pair[1];
            } else if (key == "hasRealP") {
                hasRealP = pair[1] == "true";
            } else {
                cout<<"Invalid argument: " << key << "=" << pair[1]<<endl;;
            }
        }

        void readFromFile(string argsFileName) {
            std::ifstream input( "arguments.txt" );
            for( std::string line; getline( input, line ); ) {
                vector<string> pair = parseLine(line, "=");
                setFieldFromPair(pair);
            }
            input.close();
        }

        void dumpArgs() {
            cout<<"Arguments: "<<endl;
            cout<<"Imax="<<Imax<<endl;
            cout<<"Jmax="<<Jmax<<endl;
            cout<<"rho="<<rho<<endl;
            cout<<"lineSpacing="<<lineSpacing<<endl;
            cout<<"eps="<<eps<<endl;
            cout<<"numOfAngles="<<numOfAngles<<endl;
            cout<<"numOfIterations="<<numOfIterations<<endl;
            cout<<"dx="<<dx<<endl;
            cout<<"dy="<<dy<<endl;
            cout<<"dataFile="<<dataFile<<endl;
            cout<<"numOfHeaderLines="<<numOfHeaderLines<<endl;
            cout<<"outFile="<<outFile<<endl;
            cout<<"hasRealP="<<hasRealP<<endl;
        }    
};

class Data {
    public:
        vector<vector<float>> dudt;
        vector<vector<float>> dvdt;
        vector<vector<float>> x;
        vector<vector<float>> y;
        vector<vector<float>> p;
        vector<vector<float>> pReal;
        vector<vector<float>> pNew;
        vector<vector<float>> pCount;
        float pAvgReal = 0;
        float pAvg = 0;
        float pDiff = 0;
        float pError = 0;

        Args* args;        

    public:
        Data(Args* args) {
            this->args = args;
            for(int i = 0; i < this->args->Imax; i++) {
                x.push_back(vector<float>(this->args->Jmax, 0));
                y.push_back(vector<float>(this->args->Jmax, 0));
                dudt.push_back(vector<float>(this->args->Jmax, 0));
                dvdt.push_back(vector<float>(this->args->Jmax, 0));
                p.push_back(vector<float>(this->args->Jmax, 0));
                pReal.push_back(vector<float>(this->args->Jmax, 0));
                pNew.push_back(vector<float>(this->args->Jmax, 0));
                pCount.push_back(vector<float>(this->args->Jmax, 0));
            }

            // generate some random numbers for testing, can remove it.
            srand((unsigned) time(NULL));
            for(int i = 0; i < this->args->Imax; i++) {
                for(int j = 0; j < args->Jmax; j++) {
                    x[i][j] = i*args->dx;
                    y[i][j] = j*args->dy;
                    pReal[i][j] = exp(sinf((x[i][j] + y[i][j])*1.0/(args->Imax + args->Jmax)*2*PI));
                    pAvgReal += pReal[i][j];
                }
            }
            pAvgReal = pAvgReal/args->Imax/args->Jmax;

            // calculate dpdx and dpdy, notice rho*dudt = - dpdx
            for(int i = 0; i < this->args->Imax; i++) {
                for(int j = 0; j < args->Jmax; j++) {
                    int is = (i == 0 ? 0 : i - 1);
                    int ie = (i == args->Imax - 1 ? args->Imax - 1 : i + 1);
                    int js = (j == 0 ? 0 : j - 1);
                    int je = (j == args->Jmax - 1 ? args->Jmax - 1 : j + 1);
                    dudt[i][j] = -1/args->rho*(pReal[ie][j] - pReal[is][j])/((ie - is)*args->dx) + (((float) rand()) / (float) RAND_MAX - 0.5)/10; 
                    dvdt[i][j] = -1/args->rho*(pReal[i][je] - pReal[i][js])/((je - js)*args->dy) + (((float) rand()) / (float) RAND_MAX - 0.5)/10; 
                }
            }
        }

        void readFromFile() {
            int lineNumber = 0;
            int index = 0;
            std::ifstream input( args->dataFile );
            if (!input) {
                cout<<"Could not open file: "<<args->dataFile<<endl;
            }
            for( std::string line; getline( input, line ); ) {
                cout<<line<<endl;
                lineNumber++;
                if (lineNumber > args->numOfHeaderLines) {
                    int i = index/this->args->Jmax;
                    int j = index%this->args->Jmax;
                    vector<string> pointWiseData = parseLine(line, " ");
                    readPointWiseData(pointWiseData, i, j);
                    index++;
                }
            }
            input.close();
        }

    private:
        void readPointWiseData(vector<string> pointWiseData, int i, int j) {
            x[i][j] = std::stof(pointWiseData[0]);
            y[i][j] = std::stof(pointWiseData[1]);
            dudt[i][j] = std::stof(pointWiseData[2]);
            dvdt[i][j] = std::stof(pointWiseData[3]);
            if (args->hasRealP) {
                pReal[i][j] = std::stof(pointWiseData[4]);
            }
        }


    public:
        void writePressureToFile() {
            ofstream outdata;
            outdata.open(args->outFile);
            outdata<<"TITLE=\"Pressure integrated from Omni2D\""<<endl;
            if (args->hasRealP) {
                outdata<<"VARIABLES = \"X\",\"Y\",\"Dudt\",\"DvDt\",\"Pcal\",\"Preal\",\"Count\""<<endl;
            } else {
                outdata<<"VARIABLES = "<< "\"X\",\"Y\",\"Dudt\",\"DvDt\",\"Pcal\",\"Count\""<<endl;
            }
            outdata<<"ZONE I="<<args->Imax<<" "<<"J="<<args->Jmax<<endl;

            for(int i = 0; i < args->Imax; i++) {
                for(int j = 0; j < args->Jmax; j++) {
                    if (args->hasRealP) {
                        outdata<<x[i][j]<<" "<<y[i][j]<<" "<<dudt[i][j]<<" "<<dvdt[i][j]<<" "<<p[i][j] - pAvg + pAvgReal<<" "<<pReal[i][j]<<" "<<pCount[i][j]<<endl;
                    } else {
                        outdata<<x[i][j]<<" "<<y[i][j]<<" "<<dudt[i][j]<<" "<<dvdt[i][j]<<" "<<p[i][j] - pAvg<<" "<<pCount[i][j]<<endl;
                    }
                }
            }
            outdata.close();
        }     

        void parallelLinesOdiMultiple() {
            cout<<"Processing starts"<<endl;
            for(int i = 0; i < args->numOfIterations ; i++) {
                parallellinesOdiOnce(i == args->numOfIterations - 1);
                if (args->hasRealP) {
                    cout<<"Iteration: "<<i<<" abs(p - pReal): "<<pError << " abs(pN - p):"<<pDiff<<endl;
                } else {
                    cout<<"Iteration: "<<i<< " abs(pN - p):"<<pDiff<<endl;
                }
            }
        }

        void parallellinesOdiOnce(bool isLastIeration) {
            int numOfParallelLines = max(args->Imax, args->Jmax)*1.414/args->lineSpacing;
            for(int angleId = 0; angleId < args->numOfAngles; angleId++) {
                float angle = (1.0 * angleId * 2 * PI)/ args->numOfAngles;
                float kx = cosf(angle);
                float ky = sinf(angle);
                for(int lineId = 0; lineId < numOfParallelLines; lineId++) {
                    float ic = sinf(angle)*(lineId - numOfParallelLines / 2)*args->lineSpacing + (args->Imax - 1) * 0.5;
                    float jc = cosf(angle)*(lineId - numOfParallelLines / 2)*args->lineSpacing + (args->Jmax - 1) * 0.5;
                    vector<point*> crossingPoints = getCrossingPointsOnBoundary(kx, ky, ic, jc);
                    if (crossingPoints.size() == 2) {
                        bodyIntegralAlongALine(kx, ky, ic, jc, crossingPoints[0]->x, crossingPoints[0]->y, crossingPoints[1]->x, crossingPoints[1]->y);
                    }
                }
            }
            // averaging
            pAvg = 0;
            pDiff = 0;
            pError = 0;
            for(int i = 0; i < args->Imax; i++) {
                for(int j = 0; j < args->Jmax; j++) {
                    pAvg += p[i][j];
                    if (pCount[i][j] > 0) {
                        pDiff += abs(pNew[i][j]/pCount[i][j] - p[i][j]);
                        p[i][j] = pNew[i][j]/pCount[i][j];
                        if (!isLastIeration) {
                            pNew[i][j] = 0;
                            pCount[i][j] = 0; 
                        }
                    }
                }
            }
            pAvg = pAvg/args->Imax/args->Jmax;
          
            for(int i = 0; i < args->Imax; i++) {
                for(int j = 0; j < args->Jmax; j++) {
                    pError += abs(p[i][j] - pAvg - (pReal[i][j] - pAvgReal));
                }
            }
            pDiff = pDiff/args->Imax/args->Jmax;
            pError = pError/args->Imax/args->Jmax;
        }

    public:
        vector<point*> getCrossingPointsOnBoundary(float kx, float ky, float ic, float jc) {
            // case 1, vertical to y axis
            if (abs(ky) < ZERO) {
                if(round(jc) < 0 || round(jc) > args->Jmax - 1)return vector<point*>();
                vector<point*> res{
                    new point(0, round(jc)),
                    new point(args->Imax - 1, round(jc))
                };
                if (kx < 0) {
                    reverse(res.begin(), res.end());
                }
                return res;
            }

            // case 2, vertical to x axis
            if (abs(kx) < ZERO) {
                if (round(ic) < 0 || round(ic) > args->Imax - 1)return vector<point*>();
                vector<point*> res{
                    new point(round(ic), 0),
                    new point(round(ic), args->Jmax - 1)
                };
                if (ky < 0) {
                    reverse(res.begin(), res.end());
                }
                return res;               
            }

            // case 3, inclined
            // line func: (j - jc)/ky = (i - ic)/kx
            vector<point*> res;
            // 1. crossing with line i = 0
            int j1 = round((0 - ic) / kx * ky + jc);
            if (j1 >= 0 && j1 < args->Jmax) {
                res.push_back(new point(0, j1));
            }

            // 2. crossing with line i = Imax - 1
            int j2 = round((args->Imax - 1 - ic) / kx * ky + jc);
            if (j2 >= 0 && j2 < args->Jmax) {
                res.push_back(new point(args->Imax - 1, j2));
            } 

            // 3. crossing with line j = 0
            int i1 = round((0 - jc)/ky*kx + ic);
            if (i1 >= 0 && i1 < args->Imax) {
                res.push_back(new point(i1, 0));
            }

            // 4. crossing with line j = Jmax - 1
            int i2 = round((args->Jmax - 1 - jc)/ky*kx + ic);
            if (i2 >= 0 && i2 < args->Imax) {
                res.push_back(new point(i2, args->Jmax - 1));
            }
            return res;
        }

        float bodyIntegralAlongALine(float kx, float ky, float ic, float jc, int iin, int jin, int iout, int jout) {
            // line: (i - ic)/kx = (j - jc)/ky
            if (kx*(iout - iin) + ky*(jout - jin) < 0) {
                int tmp = iout;
                iout = iin;
                iin = tmp;
                tmp = jout;
                jout = jin;
                jin = tmp;
	        }
            int ilast, jlast, inext1, jnext1, inext2, jnext2;
            ilast = iin;
            jlast = jin;

            float pint = 0;
            int choose = -1;
            pNew[ilast][jlast] += p[ilast][jlast];
            pCount[ilast][jlast] += 1;
            do {
                choose = -1;
		        if (ilast < iout) {
			        inext1 = ilast + 1; jnext1 = jlast;
		        } else if (ilast == iout) {
			        inext1 = ilast - 60000; jnext1 = jlast;
		        } else {
			        inext1 = ilast - 1; jnext1 = jlast;
		        }

                if (jlast < jout) {
			        inext2 = ilast; jnext2 = jlast + 1; 
		        } else if (jlast == jout) {
			        inext2 = ilast; jnext2 = jlast - 60000; 
		        } else {
			        inext2 = ilast; jnext2 = jlast - 1;
		        }

                // line func: (i - ic)/kx = (j - jc)/ky
                // ky*i -ic*ky - kx*j + jc*kx = 0;
                // distance= |ky*i -ic*ky - kx*j + jc*kx|/sqrt(kx^2 + ky^2)
                float d1 = abs(ky*inext1 - ky*ic - kx*jnext1 + kx*jc) / sqrt(kx*kx + ky*ky);
                float d2 = abs(ky*inext2 - ky*ic - kx*jnext2 + kx*jc) / sqrt(kx*kx + ky*ky);
                // integral marching towards minimal distance
                if (abs(d1 - d2) < 1e-5) {
                    if (inext1 >= 0 && inext1 < args->Imax) {
                        choose = 1;
                    } else if (jnext2 >= 0 && jnext2 < args->Jmax) {
                        choose = 2;
                    }
                } else {
                    if (d1 < d2 && inext1 >= 0 && inext1 < args->Imax) {
                        choose = 1;
                    } else if(d1 > d2 && jnext2 >= 0 && jnext2 < args->Jmax) {
                        choose = 2;
                    }
                }
                if (choose == 1) { 
                    pint += -args->rho*(inext1 - ilast)*args->dx*0.5*(dudt[inext1][jnext1] + dudt[ilast][jlast]);
                    pNew[inext1][jnext1] += (p[iin][jin] + pint);
                    pCount[inext1][jnext1] += 1;
                    ilast = inext1;
                } else if (choose == 2) {
                    pint += -args->rho*(jnext2 - jlast)*args->dy*0.5*(dvdt[inext2][jnext2] + dvdt[ilast][jlast]);
                    pNew[inext2][jnext2] += (p[iin][jin] + pint);
                    pCount[inext2][jnext2] += 1;
                    jlast = jnext2;                 
                }
            } while(abs(ilast - iout) + abs(jlast - jout) > 1e-1 && choose != -1);
            //if (choose == -1)cout<<"Error, wrong point"<<endl;
            return choose != -1 ? pint : 0;
        }        
};

int main(){
    Args* args = new Args();
    args->readFromFile("arguments.txt");
    args->dumpArgs();
    Data* data = new Data(args);
    data->readFromFile();
    data->parallelLinesOdiMultiple();
    data->writePressureToFile();
    delete args;
    delete data;
    return 0;
}