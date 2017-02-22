//
//  planefileformat.cpp
//  searchNoise
//
//  Created by carushi on 11/10/27.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//


#include "planefileformat.h"
namespace mousebrain
{
    const double PlaneExpression::THRESHOLD = 0.0;


    bool sort_first(const Element& first, const Element& second)
    {
        return first.first < second.first;
    }

   bool sort_distance(const pair<double, double>& first, const pair<double, double>& second)
    {
        return first.first < second.first;
    }


    bool sort_second(const Element& first, const Element& second)
    {
        return first.second < second.second;
    }

    void PlaneExpression::freeData()
    {
        for (int i = 0; i < dataGene; i++) {
            for (int j = 0; j < XNUM; j++)
                delete[] planeData[i][j];
            delete[] planeData[i];
        }
        delete[] planeData;
    }

    void PlaneExpression::clearData()
    {
        for (int i = 0; i < dataGene; i++)
            for (int j = 0; j < XNUM; j++)
                for (int k = 0; k < YNUM; k++)
                    planeData[i][j][k] = 0.0;
        return;
    }

    void PlaneExpression::reserveData()
    {
        planeData = new double**[dataGene];
        for (int i = 0; i < dataGene; i++) {
            planeData[i] = new double*[XNUM];
            for (int j = 0; j < XNUM; j++) {
                planeData[i][j] = new double[YNUM];
                for (int k = 0; k < YNUM; k++)
                    planeData[i][j][k] = 0.0;
            }
        }
        return;
    }
    void PlaneExpression::setAnnotation(int z)
    {
        class AnnotationTree tree(z);
        annotationData = tree.getAnnotationNum();
        return;
    }
    void PlaneExpression::initDataSetting(int z, bool foldchange, IDLIST* geneList)
    {
        setAnnotation(z);
        getPlaneDataFromBF(z,foldchange, geneList);
        return;
    }

    void PlaneExpression::writeBF(bool foldchange)
    {
        if (dataGene != allGene) return;
        string filename = (foldchange) ? RawFilename::numToBF(200+getz) : RawFilename::numToBF(getz);
        ofstream ofs(filename.c_str(), ios::binary);
        cout << "write" << filename << endl;
        for (int i = 0; ofs && i < allGene; i++)
            for (int j = 0; j < XNUM; j++)
                for (int k = 0; k < YNUM; k++)
                    if ( !foldchange || ( foldchange && isData(j, k) ) ) {
                        if ( foldchange ) {
                            double value = getPointData(i, j, k);
                            if (!ofs.write((char*)&value, sizeof(double)))
                                cout << "error" << endl;
                        }
                    }

        ofs.close();
        return;
    }

    void PlaneExpression::getPlaneDataFromBF(int z, bool foldchange, IDLIST* geneList)
    {
        annMax = 0;
        getz = z;
        if (geneList == NULL) {
            IDLIST tempList = IDLIST(allGene, 1);
            geneList = &tempList;
        }
        string filename = (foldchange) ? RawFilename::numToBF(200+getz) : RawFilename::numToBF(getz);
        ifstream ifs(filename.c_str(), ios::binary);
        if (!ifs) return;
        for (int i = 0, count = 0; ifs && i < allGene && count < dataGene; i++) {
            for (int j = 0; j < XNUM; j++) {
                for (int k = 0; k < YNUM; k++) {
                    if ( !foldchange || ( foldchange && isData(j, k) ) ) {
                        if ( count == 0 && isData(j, k) ) annMax++;
                        double value = 0.0;
                        ifs.read((char*)&value, sizeof(double));
                        if (allGene == dataGene || (*geneList)[i] > 0)
                            setPointData(count, j, k, value);
                    }
                }
            }
            if ((*geneList)[i] > 0) count++;
        }
        return;
    }

    void PlaneExpression::getRankVector(Data& exprAndNum, Data& rankData, bool option)
    {
        rankData.clear();
        sort(exprAndNum.begin(), exprAndNum.end(), sort_first);
        for (int i = 0; i < exprAndNum.size(); ) {
            int tail = getSameRank(i, exprAndNum);
            double truerank = (double)(tail+1+i)/2.0;
            for (; i < tail; i++) rankData.push_back(Element(truerank, exprAndNum[i].second));
        }
        exprAndNum.clear();
        if (option) sort(rankData.begin(), rankData.end(), sort_second);
        return;
    }

    /* * * * * * * * * * * * * * * * */

    void PlaneExpressionForCorrelation::getPlaneData(int z)
    {
        getz = z;
        for (int i = 0; i < (int)maxZlist.size(); i++) {
            getOneGenePlaneData(i, maxZlist[i].first, getz);
            if ( i%10000 == 0 ) cout << "a" << endl;
        }
        return;
    }

    void PlaneExpressionForCorrelation::getOneGenePlaneData(int gene, int svafilename, int z)
    {
        string str = numToStr(svafilename);
        ifstream ifs(str.c_str());
        getline(ifs,str);      getline(ifs,str);
        while ( getline(ifs, str) ) {
            Words v = myParse( str, "," );
            if ( atoi(v[2].c_str()) > z ) break;
            if ( atoi(v[2].c_str()) == z )
                setPointData( gene, atoi(v[0].c_str()), atoi(v[1].c_str()), atof(v[3].c_str()) );
        }
        ifs.close();
        return;
    }

    void PlaneExpressionForCorrelation::getMedianListNew(vector<double> &medianList)
    {
        ofstream ofs(getMVFile());
        for (int i = 0; i < dataGene; i++)
        {
            boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::median, boost::accumulators::tag::variance> > acc;
            string str = numToStr(maxZlist[i].first);
            ifstream ifs(str.c_str());
            getline(ifs,str);      getline(ifs,str);
            for (int z = -1; getline(ifs, str) && z < ZNUM; )
            {
                Words v = myParse( str, "," );
                if ( v.size() == 4 ) {
                    int tempx = atoi(v[0].c_str()), tempy = atoi(v[1].c_str()), tempz = atoi(v[2].c_str());
                    if ( tempz != z ) setAnnotation((z = tempz));
                    if ( isData(tempx, tempy) ) acc(atof(v[3].c_str()));
                }
            }
            medianList.push_back(boost::accumulators::median(acc));
            cout << i << " " << maxZlist[i].first << ":" << boost::accumulators::median(acc) << " " << boost::accumulators::variance(acc) << endl;
            ifs.close();
        }
        ofs.close();
        return;
    }
    void PlaneExpressionForCorrelation::getMedianList(vector<double> &medianList)
    {
        ifstream ifs(getMVFile());
        vector<string> v;
        fileRead(ifs, v, " :");
        while( fileRead(ifs, v, " :") ) {
            if ( v.size() == 4 ) medianList.push_back(atof(v[2].c_str()));
            else if ( v.size() > 1 ) cout << v[0] << v[1] << endl;
            else cout << "no median data error" << endl;
        }
        ifs.close();
        return;
    }

    void PlaneExpression::printDataNum() const
    {
        int count1 = 0, count2 = 0;
        for (int i = 0; i < dataGene; i++) {
            for (int x1 = 0; x1 < XNUM; x1++)
                for (int y1 = 0; y1 < YNUM; y1++) {
                    if ( getPointData(i, x1, y1) > 0.0 ) {
                    if ( isData(x1, y1) )count2++;
                        count1++;
                    }
                }
        }
        cout << "all " << count1 << " ann " << count2 << endl;
        return;
    }

    void PlaneExpressionForHGD::printDataNumOption() const
    {
        if (reduce) return;
        int count1 = 0, count2 = 0;
        for (int i = 0; i < allGene; i++) {
            for (int x1 = 0; x1 < XNUM; x1++)
                for (int y1 = 0; y1 < YNUM; y1++) {
                    if ( getPointData(i, x1, y1) > 0.0 ) {
                        if ( filenameList[i] > 0 ) count2++;
                        count1++;
                    }
                }
        }
        cout << "ann " << count1 << " havseq " << count2 << endl;
        return;
    }

    void PlaneExpressionForCorrelation::writeFoldChangeData(double gamma)
    {
        vector<double> medianList;
        getMedianList(medianList);
        cout << medianList.size() << endl;
        cout << "gene median variance" << endl;
        cout << (*medianList.begin()) << endl;
        for (int z = DATASTART; z < DATAEND; z++)
        {
            cout << z << "read start" << endl;
            setAnnotation(z);
            getPlaneData(z);
            writeBF(false);
            initDataSetting(z, false, NULL);
            for (int i = 0; i < allGene; i++)
                for (int j = 0; j < XNUM; j++)
                  for (int k = 0; k < YNUM; k++) {
                      if ( !isData(j, k) || getPointData(i, j, k) <= 0.0 )
                          setPointData(i, j, k, 0.0);
                      else {
                          setPointData(i, j, k, getPointData(i, j, k)/(double)(gamma+medianList[i]));
                          if ( getPointData(i, j, k) == 0 )
                              cout << "ketaochi.." << endl;
                      }

                  }
            writeBF(true);
        }
        return;
    }

    double PlaneExpressionForCorrelation::calcCor(Data& valueAndNumX, Data& valueAndNumY)
    {
        if ( valueAndNumX.size() != valueAndNumY.size() ) cout << "error!" << endl;
        double sum_sq_x = 0.0, sum_sq_y = 0.0, sum_coproduct = 0.0, mean_x = 0.0, mean_y = 0.0, n = 1.0;
        for (int i = 0; i < (int)valueAndNumX.size() && i < (int)valueAndNumY.size(); i++)
        {
            double x = valueAndNumX[i].first, y = valueAndNumY[i].first, sweep = (n-1.0)/n,
                delta_x = (x-mean_x), delta_y = (y-mean_y);
            sum_sq_x += delta_x * delta_x * sweep;
            sum_sq_y += delta_y * delta_y * sweep;
            sum_coproduct += delta_x * delta_y * sweep;
            mean_x += delta_x / n;
            mean_y += delta_y / n++;
        }
        if ( sum_sq_x == 0.0 || sum_sq_y == 0.0 ) return numeric_limits<double>::quiet_NaN();
        else return sum_coproduct/sqrt(sum_sq_x*sum_sq_y);
    }

    double PlaneExpressionForCorrelation::calcCorrelation(Data& data, int x1, int y1, bool rank)
    {
        Data exprAndNumX, exprAndNumY, rankDataX, rankDataY;
        for (int i = 0; i < dataGene; i++)
            if ( getPointData(i, x1, y1) > THRESHOLD && data[i].first >THRESHOLD ) {
                exprAndNumX.push_back(data[i]);
                exprAndNumY.push_back(Element(getPointData(i, x1, y1), i));
            }
        if ( exprAndNumX.size() < MINGENE ) return numeric_limits<double>::quiet_NaN();
        if (!rank) return calcCor(exprAndNumX, exprAndNumY);
        getRankVector(exprAndNumX, rankDataX, true);
        getRankVector(exprAndNumY, rankDataY, true);
        return calcCor(rankDataX, rankDataY);
    }

    double PlaneExpressionForCorrelation::calcCorrelation(int x1, int y1, int x2, int y2, bool rank)
    {
        Data exprAndNumX, exprAndNumY, rankDataX, rankDataY;
        for (int i = 0; i < dataGene; i++)
            if ( getPointData(i, x1, y1) > THRESHOLD && getPointData(i, x2, y2) > THRESHOLD ) {
                exprAndNumX.push_back(Element(getPointData(i, x1, y1), i));
                exprAndNumY.push_back(Element(getPointData(i, x2, y2), i));
            }

        if ( exprAndNumX.size() < MINGENE ) return numeric_limits<double>::quiet_NaN();
        if (!rank) return calcCor(exprAndNumX, exprAndNumY);
        getRankVector(exprAndNumX, rankDataX, true);
        getRankVector(exprAndNumY, rankDataY, true);
        return calcCor(rankDataX, rankDataY);
    }

    void PlaneExpressionForCorrelation::getDataForMaxZ()
    {
        for (int i = 0; i < (int)maxZlist.size(); i++) {
            getOneGenePlaneData(i, maxZlist[i].first, maxZlist[i].second);
            if ( i % 10000 == 0 ) cout  << "a" << endl;
        }
        return;
    }
    void PlaneExpressionForCorrelation::appCorrelation(int x1, int y1, double start, double end, vector<pair<double, double> > &dataList)
    {
        for (int i = 0; i < XNUM; i++)
            for (int j = 0; j < YNUM; j++) {
                if ( ( i == x1 && j == y1 ) || !isData(i, j) || distance(i, j, x1, y1) < start || distance(i, j, x1, y1 ) >= end ) continue;
                double data = calcCorrelation(x1, y1, i, j, true);
                if ( data != data ) data = 0.0;
                dataList.push_back(pair<double, double>(distance(x1, y1, i, j), data));
            }
        sort(dataList.begin(), dataList.end(), sort_distance);
        return;
    }
    double PlaneExpressionForCorrelation::calcMean(double correlation_min, vector<pair<double, double> > &dataList)
    {
        double mean = 0;
        for (int i = 0; i < (int)dataList.size(); i++)
        {
            mean = mean*(double)i/(double)(i+1) + dataList[i].second/(double)(i+1);
            if ( i >= 3 && mean < correlation_min ) return dataList[i].first;
        }
        return -1.0;
    }

    double PlaneExpressionForCorrelation::getHeatMapData(int x1, int y1, double correlation_min)
    {
        if ( !isData(x1, y1) ) return -1.0;
        vector<pair<double, double> > dataList;
        appCorrelation(x1, y1, 0, DIAMETER, dataList);
        double value = calcMean(correlation_min, dataList);
        if ( value < 0.0 ) {
            appCorrelation(x1, y1, DIAMETER, HEATMAX, dataList);
            value = calcMean(correlation_min, dataList);
            if ( value < 0.0 && dataList.size() != 0 )
                value = HEATMAX;
        }
        return value;
    }

    double PlaneExpressionForCorrelation::getHeatMapDataDistance(int x1, int y1, double distance_max)
    {
        double mean = 0.0; int count = 0;
        for (int i = 0; i < XNUM; i++)
            for (int j = 0; j < YNUM; j++) {
                if ( ( i == x1 && j == y1 ) || distance(x1, y1, i, j) > distance_max) continue;
                mean += calcCorrelation(x1, y1, i, j, true);
                count++;
            }
        return mean/(double)count;
    }

    void PlaneExpressionForCorrelation::setPointList(int z, bool foldchange, Table &pointList)
    {
        for (int tempz = 0; tempz < z; tempz++) {
            pointList.push_back(vector<pair<int, int> >(0));
            if ( tempz < DATASTART || tempz >= DATAEND ) continue;
            initDataSetting(tempz, foldchange, NULL);
            //printDataNum();
            for (int i = 0; i < XNUM; i++)
                for (int j = 0; j < YNUM; j++) {
                    if ( !isData(i, j) ) continue;
                    for (int k = 0, count = 0; k < dataGene; k++) {
                        if ( getPointData(k, i, j) > 0.0  ) count++;
                        if ( count >= MINGENE ) {
                            pointList[tempz].push_back(pair<int, int>(i, j));
                            break;
                        }
                    }
                }
        }
        return;
    }

    string PlaneExpressionForCorrelation::printCorData(int x1, int y1, int z1, int ann1, int x2, int y2, int z2, double correlation)
    {
        ostringstream str;
        str << x1 << " " << y1 << " " << z1 << " " << ann1 << " " << x2 << " " << y2 << " " << z2
            << " " << getAnnData(x2, y2) << " " << distance(x1, y1, z1, x2, y2, z2) << " " <<
            correlation << " " << isSameAnn(ann1, getAnnData(x2, y2));
        return str.str();
    }

    void PlaneExpressionForCorrelation::outputAllCorrelation(bool rank, bool foldchange, int z1)
    {
        string correlationfilename = RawFilename::zToFilenameFC(z1, rank, foldchange);
        ofstream ofs(correlationfilename.c_str());
        ofs << "x1 y1 z1 ann1 x2 y2 z2 ann2 distance value sameann" << endl;
        Table pointList;
        setPointList(DATASTART, foldchange, pointList);
        cout << z1 << "start!" << rank << endl;

        for (z1 = DATASTART; z1 < DATAEND; z1++) {
            cout << z1 << "reading" << endl;
            setAnnotation(z1);
            pointList.push_back(vector<pair<int, int> >(0));
            for (int x1 = XNUM-1, y1 = YNUM-1, ann1; x1 >= 0; y1--)
            {
                if ( y1 < 0 ) { x1--; y1 = YNUM-1; }
                if ( ( ann1 = getAnnData(x1, y1) ) == 0 || x1 < 0 ) continue;
                getPlaneDataFromBF(z1, foldchange, NULL);
                int count = 0;
                Data data;
                for (int i = 0; i < dataGene; i++) {
                    data.push_back(Element(getPointData(i, x1, y1), i));
                    if ( getPointData(i, x1, y1) > 0.0 ) count++;
                }
                if ( count < dataGene/100 ) continue;
                for (int z2 = z1; z2 >= DATASTART; z2--)
                    if ( pointList[z2].size() != 0 )
                    {
                        if ( z2 != z1 ) initDataSetting(z2, foldchange, NULL);
                            for (vector<pair<int, int> >::iterator it = pointList[z2].begin(); it != pointList[z2].end(); it++)
                            {
                                int x2 = it->first, y2 = it->second;
                                double correlation = calcCorrelation(data, x2, y2, rank);
                                if ( correlation == correlation ) {
                                    string output = printCorData(x1, y1, z1, ann1, x2, y2, z2, correlation);
                                    ofs << output << endl;
                                }
                            }
                    }
            cout << "(" << x1 << "," << y1  << "," << z1 << ")" << endl;
            pointList[z1].push_back(pair<int, int>(x1, y1));
            }
        }
        return;
    }

    void PlaneExpressionForCorrelation::allCalcCor(ofstream& ofs, Table& pointList, bool rank, bool foldchange)
    {
        vector<Point> randomTable;
        for (int i = 0; i < pointList.size(); i++)
            for (int j = 0; j < pointList[i].size(); j++)
                randomTable.push_back(Point(pointList[i][j].first, pointList[i][j].second, i));
        if ( randomTable.size() == 0 ) return;
        else randomCalcCor(ofs, randomTable, randomTable, rank, foldchange, 10000);
    }

    void PlaneExpressionForCorrelation::randomCalcCor(ofstream &ofs, vector<Point> &randomTable1,
                               vector<Point> &randomTable2, bool rank, bool foldchange, int count)
    {
        if ( count == 0 ) return;
        boost::mt19937 seed1((unsigned int)(time(NULL)));
        boost::uniform_real<> range1(0, randomTable1.size());
        boost::variate_generator<boost::mt19937, boost::uniform_real<> > rand1(seed1, range1);
        int x1, y1, z1, ann1, x2, y2, z2, first, second;
        first = (int)rand1();
        x1 = randomTable1[first].x, y1 = randomTable1[first].y, z1 = randomTable1[first].z;
        initDataSetting(z1, foldchange, NULL);
        ann1 = getAnnData(x1, y1);
        boost::mt19937 seed2((unsigned int)(time(NULL)));
        boost::uniform_real<> range2(0, randomTable2.size());
        boost::variate_generator<boost::mt19937, boost::uniform_real<> > rand2(seed2, range2);
        for (int tern = 0; tern < count; tern++)
        {
            Data data;
            for (int i = 0; i < dataGene; i++)
                data.push_back(Element(getPointData(i, x1, y1), i));
            second = (int)rand2(), x2 = randomTable2[second].x, y2 = randomTable2[second].y, z2 = randomTable2[second].z;
            initDataSetting(z2, foldchange, NULL);
            double correlation = calcCorrelation(data, x2, y2, rank);
            if ( correlation == correlation ) {
                string output = printCorData(x1, y1, z1, ann1, x2, y2, z2, correlation);
                ofs << output << endl;
            }
            if ( tern != count ) {
                first = (int)rand1(), x1 = randomTable1[first].x, y1 = randomTable1[first].y, z1 = randomTable1[first].z;
                initDataSetting(z1, foldchange, NULL);
                ann1 = getAnnData(x1, y1);
            }
        }
        return;
    }
    void PlaneExpressionForCorrelation::annotationCalcCor(ofstream &ofs, Table &pointList, int id, int pid, bool rank, bool foldchange)
    {
        vector<Point> target, parent, other;
        for (int i = DATASTART; i < pointList.size(); i++)
        {
            setAnnotation(i);
            for (int j = 0, ann; j < pointList[i].size(); j++)
            {
                ann = getAnnData(pointList[i][j].first, pointList[i][j].second);
                if (ann == id) target.push_back(Point(pointList[i][j].first, pointList[i][j].second, i));
                else if (ann == pid) parent.push_back(Point(pointList[i][j].first, pointList[i][j].second, i));
                else other.push_back(Point(pointList[i][j].first, pointList[i][j].second, i));
            }
        }
        int min = ( target.size() < MINPOINT ) ? (int)target.size() : MINPOINT;
        randomCalcCor(ofs, target, target, rank, foldchange, min);
        randomCalcCor(ofs, target, parent, rank, foldchange, (parent.size() < min) ? (int)parent.size():min);
        randomCalcCor(ofs, target, other, rank, foldchange, min);
        return;
    }

    void PlaneExpressionForCorrelation::outputCorWithRandomization(bool rank, bool foldchange)
    {
        const char* filename = getCorFile(foldchange, rank);
        ofstream ofs(filename, ios::app);
        cout << "all " << filename << endl;
        //ofs << "x1 y1 z1 ann1 x2 y2 z2 ann2 distance value sameann" << endl;
        Table pointList;
        setPointList(DATAEND, foldchange, pointList);
        if ( pointList.size() < DATAEND ) { cout << "data0" << endl; return; }
        //allCalcCor(ofs, pointList, rank, foldchange);
        ofs.close();
        ifstream ifs(RawFilename::structureFile.c_str());
        Words v;
        AnnotationTree annotation(0);
        annotation.setTree();
        fileRead(ifs, v, ",");
        //ofs.open("/Users/username/program/atlas/data/correlation/annall.txt");
        while (fileRead(ifs, v, ",")) {
            cout << v.size() << endl;
            if ( v.size() < 8 || atoi(v[6].c_str()) < 37 || atoi(v[6].c_str()) >= 40 ) continue;
            string annfile = RawFilename::zToFilename(atoi(v[6].c_str()), rank);
            int parentId = annotation.getParentId(v[2]);
                if ( parentId >= 0 ) {
                ofs.open(annfile.c_str());
                cout << v[1] << " " << v[6] << " " << v[2] << parentId << endl;
                ofs << "% " << v[1] << " " << v[6] << " " << v[2] << parentId << endl;
                ofs << "x1 y1 z1 ann1 x2 y2 z2 ann2 distance value sameann" << endl;
                annotationCalcCor(ofs, pointList, atoi(v[6].c_str()), parentId, rank, foldchange);
                ofs.close();
            }
        }
        annotation.freeTree();
        return;
    }

    void PlaneExpressionForCorrelation::outputCorForClustering()
    {
        const char* filename = RawFilename::clusterFile.c_str();
        ofstream ofs(filename, ios::trunc);
        cout << "all " << filename << endl;
        //ofs << "x1 y1 z1 ann1 x2 y2 z2 ann2 distance value sameann" << endl;
        Table pointListbefore, pointList;
        setPointList(DATAEND, true, pointListbefore);

        if ( pointListbefore.size() < DATAEND ) { cout << "data0" << endl; return; }
        int count = 0;
        for (int i = 0; i < pointListbefore.size(); i++) {
            pointList.push_back(vector<pair<int, int> >(0));
            for (vector<pair<int, int > >::iterator it = pointListbefore[i].begin(); pointListbefore[i].size() != 0 && it != pointListbefore[i].end(); it++, count++) {
                if ( count%100 == 0 ) pointList[i].push_back(pair<int,int>(it->first, it->second));
            }
        }
        cout << "eraseend" << endl;
        for (int i = 0; i < pointList.size(); i++) {
            if ( pointList[i].size() == 0 ) continue;
            setAnnotation(i);
            for (vector<pair<int, int> >::iterator it = pointList[i].begin(); it != pointList[i].end(); it++) {
                ofs << ",( " << it->first << " " << it->second << " " << i << " " << getAnnData(it->first, it->second) << " ) ";
                cout << "( " << it->first << " " << it->second << " " << i << " " << getAnnData(it->first, it->second) << " ) " << endl;
            }
        }
        ofs << endl;
        cout << endl;

        for (int z1 = DATASTART, count = 0; z1 < DATAEND; z1++) {
            cout << z1 << "reading" << endl;
            if ( pointList[z1].size() == 0 ) continue;
            for (vector<pair<int, int> >::iterator it = pointList[z1].begin(); it != pointList[z1].end(); it++) {
                initDataSetting(z1, true, NULL);
                int x1 = it->first, y1 = it->second, ann1 = getAnnData(x1, y1);
                ofs << "( " << x1 << " " << y1 << " " << z1 << " " << ann1 << " ),";
                cout << "( " << x1 << " " << y1 << " " << z1 << " " << ann1 << " ),";
                Data data;
                for (int i = 0; i < dataGene; i++)
                    data.push_back(Element(getPointData(i, x1, y1), i));
                for (int z2 = DATASTART; z2 <= z1; z2++)
                    if ( pointList[z2].size() != 0 )
                    {
                        if ( z1 != z2 ) initDataSetting(z2, true, NULL);
                        for (vector<pair<int, int> >::iterator it2 = pointList[z2].begin(); it2 != pointList[z2].end(); it2++) {
                            int x2 = it2->first, y2 = it2->second;
                            if ( z2 == z1 && x2 == x1 && y2 == y1 ) break;
                            double correlation = calcCorrelation(data, x2, y2, true);
                            if ( correlation == correlation )
                                ofs << 1.0-correlation << ",";
                            else
                                ofs << 2.0 << ",";
                        }
                    }
                for (int i = 0; i < 255-count; i++) {
                    if ( i == 255-count-1 ) ofs << "NA" << endl;
                    else ofs << "NA,";
                }
                count++;
            }
        }
        return;
    }

    void PlaneExpressionForCorrelation::outputCorWithAnnotate(bool rank, int z, double gamma)
    {
        const char* filename = getCorFile(true, rank);
        string str = ( z >= 0 ) ? RawFilename::zgammaToFilename(z, gamma, rank): filename;
        ofstream ofs(str.c_str(), ios::out);
        cout << str << endl;

        ofs << "\"( x1 y1 ann )\" \"( x2 y2 ann )\" distance value samaann" << endl;
        for( int i = XNUM-1; i >= 0; i-- )
            for ( int j = YNUM-1; j >= 0; j-- )
                for ( int k = i; k >= 0; k--)
                {
                    if ( !isData(i, j) ) continue;
                    int l = ( k == i ) ? j-1 : XNUM-1;
                    for ( ; l >= 0; l-- ) {
                        if ( !isData(k, l) ) continue;
                        double correlation = calcCorrelation(i, j, k, l, rank);
                        if ( correlation == correlation )
                            ofs << "\"( " << i << " " << j << " " << getAnnData(i, j) << " )\" \"( " << k << " " << l << " " << getAnnData(k, l) << ")\" " << distance(i, j, k, l) << " " << correlation << " " << isSameAnn(i, j, k, l) << endl;
                    }
                }
        ofs.close();
        return;
    }


    void PlaneExpression::printData() const
    {
        for (int j = 0; j < XNUM; j++)
            for (int k = 0; k < YNUM; k++) {
                cout << j << " " << k << ":" << endl;
                for (int i = 0; i < dataGene; i++)
                    cout << getPointData(i, j, k) << " ";
                cout << endl;
            }
        return;
    }

    void PlaneExpression::setFoldChangeAgain(double gamma, double ***original)
    {
        for (int i = 0; i < dataGene; i++)
        {
            boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::median> > acc;
            for (int j= 0; j < XNUM; j++)
                for (int k = 0; k < YNUM; k++)
                    if ( isData(j,k) ) acc(original[i][j][k]);
            double tmpMedian = boost::accumulators::median(acc);
            for (int j = 0; j < XNUM; j++)
                for (int k = 0; k < YNUM; k++)
                    setPointData(i, j, k, (original[i][j][k]/(tmpMedian+gamma)));
        }
        return;
    }

    void PlaneExpression::getMedianVariance(double ***original)
    {
        vector<pair<double, double> > list;
        for (int i = 0; i < dataGene; i++)
        {
            boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::median, boost::accumulators::tag::variance> > acc;
            for (int j= 0; j < XNUM; j++)
                for (int k = 0; k < YNUM; k++)
                    if ( isData(j,k) ) acc(original[i][j][k]);
            double tmpMedian = boost::accumulators::median(acc), tmpVariance = boost::accumulators::variance(acc);
            list.push_back(pair<double, double> (tmpMedian, tmpVariance));
        }
        string str = RawFilename::medianAndVariance(getz);
        ofstream ofs(str.c_str());
        ofs << "gene: median variance" << endl;
        for (int i = 0; i < (int)list.size(); i++)
            ofs << i << ": " << list[i].first << " " << list[i].second << endl;
        ofs.close();
        return;
    }

    void PlaneExpressionForCorrelation::heatAndCorrelate(double threshold)
    {
        double ***original = getAllPlaneData();
        reserveData();
        setAnnotation(getz);
        cout << "getAnnotation!" << endl;
        getMedianVariance(original);
        for (double gamma = START; gamma < END; gamma += STEP)
        {
            cout << "gamma" << gamma << endl;
            setFoldChangeAgain(gamma, original);
            outputCorWithAnnotate(true, getz, gamma);
            outputCorWithAnnotate(false, getz, gamma);
        }
        freeData();
        setAllPlaneData(original);
    }

    /* * * * * * * * * * * * * */

    const int PlaneExpressionForHGD::BIN = 10;

    int PlaneExpressionForHGD::countData()
    {
        int count = 0;
        for (IDLIST::iterator it = filenameList.begin(); it != filenameList.end(); it++)
            if (*it > 0) count++;
        return count;
    }
    void PlaneExpressionForHGD::initRes(int z, bool foldchange)
    {
        readFilenameHavSeq();
        dataGene = allGene;
        reserveData();
        initDataSetting(z, foldchange, &filenameList);
        reduce = false;
    }

    void PlaneExpressionForHGD::initResSaveMem(int z, bool foldchange)
    {
        readFilenameHavSeq();
        dataGene = countData();
        reserveData();
        initDataSetting(z, foldchange, &filenameList);
        reduce = true;
        return;
    }

    double PlaneExpressionForHGD::hyperGeographicDistribution(Data &expr, int bin, bool revision, bool upper)
    {
        int N = (int)expr.size(), n = (int)expr.size()/bin, m, k;
        getM(m, k, bin, expr, revision, upper);
        if ( m > N ) m = N;
        return sumProbability(N, n, m, k, revision);
    }

    int PlaneExpressionForHGD::threeBaseNum(string str)
    {
        int num = 0;
        for (int i = 0; i < str.length(); i++)
            switch(i) {
                case 0: num += (encodeChar(str[i])-1); break;
                case 1: num += (encodeChar(str[i])-1)*4; break;
                case 2: num += (encodeChar(str[i])-1)*16; break;
                default: i = (int)str.length();
            }
        return num;
    }

    void PlaneExpressionForHGD::getCount(double *allArray, double *upperArray, SeedEffect &binGene)
    {
        ifstream ifs(RawFilename::threebase.c_str(), ios::in | ios::binary);
        for (int i = 0, isid, count; i < dataGene; i++)
        {
            ifs.read((char*)&isid, sizeof(int));
            SeedEffect::iterator it = binGene.find(isid);
            if (RawFilename::mask) {
            double array[THREEBASE], sum;
            for (int j = 0; j < THREEBASE; j++)
            {
                ifs.read((char*)&count, sizeof(int));
                array[j] = count;
                sum += (double)count;
            }
            if ( it != binGene.end() ) {
                for (int j = 0; j < THREEBASE; j++) {
                    if ( it->second == UPPER )
                        upperArray[j] += array[j]/sum;
                    allArray[j] += array[j]/sum;
                }
            }
            } else {
             for (int j = 0; j < THREEBASE; j++)
             {
                 ifs.read((char*)&count, sizeof(int));
                 if ( it != binGene.end() ) {
                     if ( it->second == UPPER )
                         upperArray[j] += (double)count;
                     allArray[j] += (double)count;
                 }
             }
            }

        }
        return;
    }

    double PlaneExpressionForHGD::getCoefficientOfM(SeedEffect &binGene, string &seed, bool revision)
    {
        if ( seed.length() != GeneSequence::SEEDLENGTH || !revision ) return 1.0;
        ofstream ofs(getDBFile(), ios::app);
        double *allArray = new double[THREEBASE], *upperArray = new double[THREEBASE],
                    coefficient = 1.0, numerator, denominator;
        fill(allArray, allArray+THREEBASE, 0.0);
        fill(upperArray, upperArray+THREEBASE, 0.0);
        getCount(allArray, upperArray, binGene);
        mutexLock();
        for (int i = (int)seed.length()-3, num; i >= 0; i--) /* N(AAA) */
        {
            num = threeBaseNum(seed.substr(i,3));
            coefficient *= (upperArray[num]+0.1)/(allArray[num]+0.1);
            ofs << "N(AAA) (" << upperArray[num] << "+0.1)/(" << allArray[num] << "+0.1)" << endl;
        }
        for (int i = (int)seed.length()-3, num; i >= 0; i--)
        {
            numerator = denominator = 0.1;
            if ( i == (int)seed.length()-3 )
                for (int j = 0; j < THREEBASE; numerator += allArray[j], denominator += upperArray[j], j++) ; /* N(NNN) */
            else {
                num = threeBaseNum("A"+seed.substr(i+1, 2));
                for (int j = 0; j < 4; numerator += allArray[num+j], denominator += upperArray[num+j], j++) ; /* N(NAA) */
            }

            if ( i == (int)seed.length()-3 )
                ofs << "N(NNN) (" << numerator << ")/(" << denominator << ")" << endl;
            else
                ofs << "N(N" << seed.substr(i+1,2) << ") (" << numerator << ")/(" << denominator << ")"<< endl;
            coefficient *= numerator/denominator;
        }
        if ( coefficient != coefficient ) coefficient = 1.0;
        delete[] upperArray;
        delete[] allArray;
        ofs << "coefficient" << coefficient << endl;
        mutexUnlock();
        return coefficient;
    }

    void PlaneExpressionForHGD::getM(int &m, int &k, int bin, Data &expr, bool revision, bool upper)
    {
        int i;
        m = k = 0;
        SeedEffect binGene;
        if ( upper ) {
            for (i = 0; i < (int)expr.size()-(int)expr.size()/bin; i++) {
                if ( seedEffect[filenameList[expr[i].second]] == 0 ) { m++; }
                binGene.insert(pair<int, int>(expr[i].second, LOWER));
            }
            for (; i < (int)expr.size(); i++) {
                if ( seedEffect[filenameList[expr[i].second]] == 0 ) { m++; k++; }
                binGene.insert(pair<int, int>(expr[i].second, UPPER));
            }
        } else {
            for (i = 0; i < (int)expr.size()/bin; i++) {
                if ( seedEffect[filenameList[expr[i].second]] != 0 ) { m++; k++; }
                binGene.insert(pair<int, int>(expr[i].second, UPPER));
            }
            for (; i < (int)expr.size(); i++) {
                if ( seedEffect[filenameList[expr[i].second]] != 0 ) { m++; }
                binGene.insert(pair<int, int>(expr[i].second, LOWER));
            }
        }

        m = (int)((double)m*getCoefficientOfM(binGene, geneSequence->miRNATarget, revision));
        if ( m < k ) m = k;
        return;
    }

    void PlaneExpressionForHGD::readFilenameHavSeq()
    {
        ifstream ifs(RawFilename::filename7.c_str());
        string str;
        int count = 0, count1 = 0;
        IDLIST newfilenameList;
        while ( getline(ifs, str) )
            newfilenameList.push_back(atoi(str.c_str()));
        for (IDLIST::iterator it = filenameList.begin(); it != filenameList.end(); it++)
            if ( find(newfilenameList.begin(), newfilenameList.end(), (*it)) == newfilenameList.end() )
            { (*it) = -1; count++; }
        for (IDLIST::iterator it = filenameList.begin(); it != filenameList.end(); it++)
            if ( *it > 0 ) count1++;
        cout << "-target datasize" << allGene << "=" << filenameList.size() << "= " << count << " (no seq) +" << count1 << " (haveseq)" << endl;
        return;
    }

    int PlaneExpressionForHGD::maxSeedGene(Data& expr, bool upper)
    {
        int m = 0;
        for (int i = 0; i < (int)expr.size(); i++)
            if ( seedEffect[filenameList[expr[i].second]] == 0 ) m++;
        if ( !upper ) m = (int)expr.size()-m;
        return m;
    }

    double PlaneExpressionForHGD::getMinPvalue(int x1, int y1, bool upper)
    {
        double min = 1.0;
        if ( !isData(x1, y1) || !reduce) return min;
        Data expr;
        for (int i = 0, count = 0; i < dataGene; i++) {
            for (; count < allGene && filenameList[count] <= 0; count++) ;
            if (getPointData(i, x1, y1) > THRESHOLD)
                expr.push_back(Element(getPointData(i, x1, y1), count));
            count++;
        }
        sort(expr.begin(), expr.end(), sort_first);
        if ( expr.size()/BIN >= MINGENE2 ) {
            int N = (int)expr.size(), m = maxSeedGene(expr, upper), k = 0;
            //cout << "N m : " << N << " " << m << endl;
            for (int n = 0; n < N/2; n++) {
                if ( upper && seedEffect[filenameList[expr[N-1-n].second]] == 0 ) k++;
                else if ( !upper && seedEffect[filenameList[expr[n].second]] != 0 ) k++;
                if ( (n+1)%(N/BIN) == 0 && n >= MINGENE) {
                    //cout << "n k : " << n << " " << k << endl;
                    double tempp = sumProbability(N, n, m, k, false);
                    if (tempp < min) min = tempp;
                }
            }
        }
        return min;
    }

    vector<double> PlaneExpressionForHGD::outputPvalue(int x, int y)
    {
        vector<double> data;
        if (!reduce) return data;
        ofstream ofs(RawFilename::pvalueFile.c_str(), ios::app);
        if ( !isData(x, y) ) { ofs << endl; return data; }
        Data expr;
        for (int i = 0, count = 0; i < dataGene; i++) {
            for (; count < allGene && filenameList[count] <= 0; count++) ;
            if (getPointData(i, x, y) > THRESHOLD)
                expr.push_back(Element(getPointData(i, x, y), count));
            count++;
        }
        sort(expr.begin(), expr.end(), sort_first);
        if ( expr.size()/20 >= MINGENE2 ) {
            for (int i = BIN; i >= 2; i--)
                data.push_back(hyperGeographicDistribution(expr, i, true, false));
            for (int i = BIN; i >= 2; i--)
                data.push_back(hyperGeographicDistribution(expr, i, false, false));
        }
        return data;
    }

}
