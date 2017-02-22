//
//  psoutput.cpp
//  searchNoise
//
//  Created by carushi on 11/09/17.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_psoutput_cpp
#define searchNoise_psoutput_cpp

#include "filename.h"

namespace mousebrain
{
    void Paint::setRGB( double *r, double *g, double *b, double data, double max )
    {
        if ( data < 0.0 ) { *r = *g = *b = 0.8; return; }
        *r = *g = *b = 0.0;
        if ( max <= 1.0 )
            *r = *g = *b = data;
        else {
            *r = data/max;
            *g = ( data > max/2.0 ) ? data/(max/2.0)-1.0 : data/(max/2.0);
            *b = ( data > max/3.0 ) ? ( ( data > max/3.0 *2.0 ) ? data/max*3.0-2.0 : data/max*3.0-1.0 ) : data/max*3.0;
        }
        return;
    }

    void Paint::setRGBmono( double *r, double *g, double *b, double data, double max )
    {
        *r = *g = *b = 0.0;
        if ( data != data || data <= 0.0 ) return;
        *r = *g = *b = data/(1.0+max);
        return;
    }

    void Paint::setPvalue( double *r, double *g, double *b, double data, double bon, int ann )
    {
        *r = *g = *b = 0.0;
        if ( ann != 0 ) {
            *r = *g = *b = 1.0;
            if ( data == data && data >= 0.0 ) {
                if ( data > 1.0 ) *r = *g = *b = 0.5;
                else if ( bon > 1.0 )  { if ( data*bon < Filename::PVALUE ) *r = *g = 0.0; *b = 1.0; }
                else { *r = data; *g = 1.0-sqrt((0.5-data)*(0.5-data)); *b = 1.0-data; }
            }
        }
        return;
    }

    void Paint::setRGBheatmap( double *r, double *g, double *b, double data, double max)
    {
        if ( data < 1.0 ) *r = *g = *b = 0.0;
        else if ( data >= 1.0 ){
            if ( data > max ) data = max;
            *r = data/max;
            *b = 1.0-sqrt((0.5-data/max)*(0.5-data/max));
            *g = 1.0-data/max;
        }
        return;
    }

    void Paint::psoutput(const char* filename, class SmoothExpression &filedata)
    {
        string outputfile = RawFilename::psFilename(filename, -1);
        ofstream ofs(outputfile.c_str());
        ofs << "%!\ngsave" << endl;

        for (int h = 0; h < Filename::ZNUM; h++) {
            for (int i = 0; i < Filename::XNUM; i++) {
                for (int j = 0; j < Filename::YNUM; j++) {
                    ofs << 20 + i + (h%6)*100 << " " << 750 - j -(h/6)*50<< " moveto" << endl;
                    double r = 0.0, g = 0.0, b = 0.0;
                    if ( strcmp(filename, "AtlasAnnotation200") == 0 )
                        setRGB( &r, &g, &b, filedata.getData(i, j, h), 232 );
                    else
                        setRGB( &r, &g, &b, filedata.getData(i, j, h), filedata.getMax());
                    ofs << r << " " << g << " " << b << " setrgbcolor" << endl;
                    ofs << 20 + i + 1+(h%6)*100 << " " << 750-j -(h/6)*50 << " lineto" << endl;
                    ofs << "stroke" << endl;
                }
            }
        }
        ofs << "showpage\ngrestore" << endl;
        ofs.close();
        return;
    }

    void Paint::psoutput(const char* filename, class PlaneExpressionForCorrelation &filedata, double correlation_min)
    {
        double r = 0.5, g = 0.5, b = 0.5, data = 0.0;
        //string outputfile = psFilename(heatmapAndGamma(filedata.getz, (0.2)).c_str(), correlation_min);
        string outputfile = RawFilename::psFilename("heatmapwithAllDataRGB", correlation_min);
        ofstream ofs(outputfile.c_str());
        if ( !ofs ) { cout << outputfile << "error" << endl; return; }
        ofs << "%!\ngsave" << endl;

        for (int h = 0; h < Filename::ZNUM; h++) {
            filedata.initDataSetting(h, true, NULL);
            for (int i = 0; i < Filename::XNUM; i++) {
                for (int j = 0; j < Filename::YNUM; j++) {
                    ofs << 20 + i+(h%6)*100 << " " << 750-j-(h/6)*50 << " moveto" << endl;
                    if ( filedata.isData(i, j) ) {
                        data = filedata.getHeatMapData(i, j, correlation_min);
                        setRGBheatmap( &r, &g, &b, data, HEATMAX);
                    } else {
                        r = g = b = 0.5;
                    }
                    if ( data > 0.0 )
                        cout << i << " " << j << " " << h << " " << data << endl;
                    ofs << r << " " << g << " " << b << " setrgbcolor" << endl;
                    ofs << 20 + i + 1+(h%6)*100 << " " << 750-j -(h/6)*50 << " lineto" << endl;
                    ofs << "stroke" << endl;
                }
            }
        }
        for (int chart = 0; chart < 10; chart++) {
            setRGBheatmap( &r, &g, &b, (double)HEATMAX/10.0*(double)(chart+1), HEATMAX);
            drawColorChart(ofs, r, g, b, chart);
        }

        ofs << "showpage\ngrestore" << endl;
        ofs.close();
        return;
    }
    void Paint::drawPs(ofstream &ofs, int x1, int y1, int z1, double r, double g, double b)
    {
        ofs << 20 + x1 + (z1%6)*100 << " " << 750 - y1 -(z1/6)*50<< " moveto" << endl;
        ofs << r << " " << g << " " << b << " setrgbcolor" << endl;
        ofs << 20 + x1 + 1+(z1%6)*100 << " " << 750- y1 -(z1/6)*50 << " lineto" << endl;
        ofs << "stroke" << endl;
        return;
    }
    void Paint::drawColorChart(ofstream &ofs, double r, double g, double b, int num)
    {
        ofs << 20+num+(DATAEND%6)*100 << " " << 300  << " moveto" << endl;
        ofs << r << " " << g << " " << b << " setrgbcolor" << endl;
        ofs << 20+num+(DATAEND%6)*100 << " " << 295 << " lineto" << endl;
        ofs << "stroke" << endl;

    }
    double Paint::getValue(int type, int x, int y, class PlaneExpressionForHGD &filedata)
    {
        switch (type) {
            case 0: case 1:
                return filedata.getMinPvalue(x, y, true);
            case 2: { return filedata.getMinPvalue(x, y, true); }
            case 3: { return filedata.getMinPvalue(x, y, false); }
            case 4: { return filedata.getMinPvalue(x, y, true); }
        }
        return 0.0;
    }
    void Paint::setColor(int type, double &r, double &g, double &b, double bon, double data, int ann)
    {
        switch (type) {
            case 0: case 1: {
                double tempbon = ( type ) ? 1.0 : bon;
                setPvalue( &r, &g, &b, data, tempbon, ann );
                return;
            }
            default:
                setPvalue( &r, &g, &b, data, 1.0, ann);
                return;
        }
    }
    void Paint::psoutput(string filename, class PlaneExpressionForHGD &filedata, int miRNANo, int bon)
    {
        double r = 0.0, g = 0.0, b = 0.0;
        filedata.getz = 0;
        vector<string>filenameList;
        for (int type = 0; type < 5; type++)
            filenameList.push_back(RawFilename::psFilename(filename.c_str(), ((double)filedata.getz+100*type)));
        for (int h = DATASTART; h < DATAEND; h++) {
            filedata.initDataSetting(h, true, NULL);
            cout << h << " " << filedata.annMax << endl;
            for (int type = 0; type < 5; type++) {
                ofstream ofs;
                if ( h == DATASTART ) {
                    ofs.open(filenameList[type].c_str(), ios::trunc);
                    ofs << "%!\ngsave" << endl;
                }
                else
                    ofs.open(filenameList[type].c_str(), ios::app);
                for (int i = 0; i < Filename::XNUM; i++)
                    for (int j = 0; j < Filename::YNUM; j++) {
                        setColor(type, r, g, b, bon, getValue(type, i, j, filedata), filedata.getAnnData(i, j) );
                        drawPs(ofs, i, j, h, r, g, b);
                    }
                ofs.close();
            }
        }
        for (int type = 0; type < 5; type++) {
            ofstream ofs(filenameList[type].c_str(), ios::app);
            for (int chart = 0; chart <= 10; chart++) {
                setColor(type, r, g, b, bon, 1.0/10.0*((double)chart+1.0), 1);
                drawColorChart(ofs, r, g, b, chart);
            }
            ofs << "showpage\ngrestore" << endl;
        }
        return;
    }

}
#endif
