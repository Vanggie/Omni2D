#include <fstream>
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
        }    
};

class Data {
    public:
        vector<vector<float>> dudt;
        vector<vector<float>> dvdt;
        vector<vector<float>> x;
        vector<vector<float>> y;
        vector<vector<float>> p;
        vector<vector<float>> pNew;
        vector<vector<float>> pCount;

        Args* args;        

    public:
        Data(Args* args) {
            this->args = args;
            for(int i = 0; i < this->args->Imax; i++) {
                x.push_back(vector<float>(this->args->Jmax));
                y.push_back(vector<float>(this->args->Jmax));
                dudt.push_back(vector<float>(this->args->Jmax));
                dvdt.push_back(vector<float>(this->args->Jmax));
                p.push_back(vector<float>(this->args->Jmax));
                pNew.push_back(vector<float>(this->args->Jmax));
                pCount.push_back(vector<float>(this->args->Jmax));
            }

            // generate some random numbers for testing, can remove it.
            srand((unsigned) time(NULL));
            for(int i = 0; i < this->args->Imax; i++) {
                for(int j = 0; j < args->Jmax; j++) {
                    dudt[i][j] = rand()/(float (rand()));
                    dvdt[i][j] = rand()/(float (rand()));
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
        }


    public:
        void writePressureToFile() {
            ofstream outdata;
            outdata.open(args->outFile);
            for(int i = 0; i < args->Imax; i++) {
                for(int j = 0; j < args->Jmax; j++) {
                    outdata<<x[i][j]<<" "<<y[i][j]<<" "<<dudt[i][j]<<" "<<dvdt[i][j]<<" "<<p[i][j]<<" "<<pCount[i][j]<<endl;
                }
            }
            outdata.close();
        }     

        void parallelLinesOdiMultiple() {
            cout<<"Processing starts"<<endl;
            for(int i = 0; i < args->numOfIterations ; i++) {
                parallellinesOdiOnce(i == args->numOfIterations -1);
                cout<<"Iteration: "<<i<<endl;
            }
        }

        void parallellinesOdiOnce(bool isLastIeration) {
            int numOfParallelLines = max(args->Imax, args->Jmax)*1.414/args->lineSpacing;
            for(int angleId = 0; angleId < args->numOfAngles; angleId++) {
                float angle = (angleId * 2 * PI)/ args->numOfAngles;
                float kx = cosf(angle);
                float ky = sinf(angle);
                for(int lineId = 0; lineId < numOfParallelLines; lineId++) {
                    float ic = sinf(angle)*(lineId - numOfParallelLines / 2)*args->lineSpacing + (args->Imax - 1) * 0.5;
                    float jc = cosf(angle)*(lineId - numOfParallelLines / 2)*args->lineSpacing + (args->Jmax - 1) * 0.5;
                    vector<point*> crossingPoints = getCrossingPointsOnBoundary(kx, ky, ic, jc);
                    if (crossingPoints.size() >= 2) {
                        bodyIntegralAlongALine(kx, ky, ic, jc, crossingPoints[0]->x, crossingPoints[0]->y, crossingPoints[1]->x, crossingPoints[1]->y);
                    }
                }
            }
            // averaging
            for(int i = 0; i < args->Imax; i++) {
                for(int j = 0; j < args->Jmax; j++) {
                    if (pCount[i][j] > 0) {
                        p[i][j] = pNew[i][j]/pCount[i][j];
                        if (!isLastIeration) {
                            pNew[i][j] = p[i][j];
                            pCount[i][j] = 1; 
                        }
                    }
                }
            }
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
                swap(iout, iin);
                swap(jout, jin);
	        }
            int ilast, jlast, inext1, jnext1, inext2, jnext2;
            ilast = iin;
            jlast = jin;

            float pint = 0;
            bool flag = 0;
            do {
                flag = 0;
		        if (ilast < iout) {
			        inext1 = ilast + 1; jnext1 = jlast;
		        } else if (ilast == iout) {
			        inext1 = ilast - 60000; jnext1 = jlast;
		        } else {
			        inext1 = ilast - 1; jnext1 = jlast;
		        }

                if (jlast<jout) {
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
                if (d1 <= d2 && inext1 >= 0 && inext1 < args->Imax) { 
                    pint += -args->rho*(inext1 - ilast)*args->dx*0.5*(dudt[inext1][jnext1] + dudt[ilast][jlast]);
                    pNew[inext1][jnext1] += p[iin][jin] + pint;
                    pCount[inext1][jnext1]++;
                    ilast = inext1;
                    flag = 1;
                } else if (d1 > d2 && jnext2 >= 0 && jnext2 < args->Jmax) {
                    pint += -args->rho*(jnext2 - jlast)*args->dy*0.5*(dvdt[inext2][jnext2] + dvdt[ilast][jlast]);
                    pNew[inext2][jnext2] += p[iin][jin] + pint;
                    pCount[inext2][jnext2]++;
                    jlast = jnext2;
                    flag = 1;                    
                }

            } while(abs(ilast - iout) + abs(jlast - jout) > 1e-5 && flag);
            return flag ? pint : 0;
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