/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2005 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
 * Tobias Ringström, and Jonas Samuelsson.                                   *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*!
  \file
  \brief Plotting of data in PDF format
  \author Tony Ottosson

  $Revision$

  $Date$ 
*/

#ifndef __plot_pdf_h
#define __plot_pdf_h

#ifdef HAVE_HARU

#include "itpp/graphics/libharu.h"

#include "itpp/base/vec.h"
#include "itpp/base/mat.h"
#include "itpp/base/array.h"
#include "itpp/base/converters.h"


/*
TODO:


* Logarithmic scale (x and/or y)
* legend
* line stiles (dotted, dash-dotted)
* markers at data-points
* line width
* different fonts
* save plot and variables in it-file format to be able to post-process plot.
*/



namespace itpp {


  enum LineStyle {SOLID, DASHED, DOTTED};

  enum MarkerType {NONE, SQUARE};

  /*!
    \brief Plotting data in PDF output
    \author Tony Ottosson
  
    This is an experimental feature class that will change in the future. Please
    use at your own risk.

    Examples:
    \code
    #include "itpp/itbase.h"
  
    int main (){
    PlotPDF plot("plot.pdf");
  
    vec x=linspace(0,1.0, 100);
    vec y= sin(2*pi*x);
  
    plot.set_plot(x,y);
    plot.set_axis(0, 1.0, -1.0, 1.0);
  
    plot.save_pdf();
    }
    \endcode
  */
  class PlotPDF {
  public:
    //! Constructor with filename input
    PlotPDF(const std::string &filename);
  
    //! Save the plot in PDF
    void save_pdf();

    //! Save plot in generic binary format
    void save();
    //! Load plot (generic binary format)
    void load();


    //! Input a plot object (x,y) in vectors. returns plot number
    int set_plot(const vec &x, const vec &y);
    //! Set line style
    void set_line_style(const int plot_number, const LineStyle style);
    //! Set line marker
    void set_line_marker(const int plot_number, const MarkerType marker);
    //! set line width
    void set_line_width(const int plot_number, const double width);

    //! Set range of axis
    void set_axis(const double xmin, const double xmax, const double ymin, const double ymax);
    //! Set grid
    void grid_on() { grid = true; }
    //! Turn off grid
    void grid_off() { grid = false; }
    
    //! Set X-label
    void set_xlabel(const std::string &text) { XLabel = text; }
    //! Set Y-label
    void set_ylabel(const std::string &text) { YLabel = text; }
    //! Set title
    void set_title(const std::string &text){ Title = text; }

  private:

    // Convert cm to pts
    double cm2pt(const double cm) { return (cm/2.54)*72;}

    // convert (x,y) coordinates to absolute values given in pts
    double x2pt(const double x) { return (x_origo + (x-XMin)*dx); }
    double y2pt(const double y) { return (y_origo + (y-YMin)*dy); }

    // convert absolute pt coordintates to (x,y) coordinates
    double pt2x(const double pt) { return ((pt-x_origo)/dx + XMin); }
    double pt2y(const double pt) { return ((pt-y_origo)/dy + YMin); }
  
    // returns true if coordinate is within current axis range
    bool in_range(const double x, const double y);

    // Move, line_to and stroke 
    void draw_line(const double x0, const double y0, const double x1, const double y1);

    void move_to(const double x, const double y);
    void line_to(const double x, const double y);
    void stroke_line();
    void draw_line_style(const LineStyle style);
    void draw_line_width(const double width);

    // draw a marker at position (x,y)
    void draw_marker(const double x, const double y, const MarkerType marker);

    void set_font(const std::string &fontname, int fontsize);

    void draw_text(const double x, const double y, const std::string &text);
    void draw_text_centered(const double x, const double y, const std::string &text);
    void draw_text_right_justified(const double x, const double y, const std::string &text);

    void draw_vertical_text_centered(const double x, const double y, const std::string &text);

    void draw_xticks();
    void draw_yticks();
    void draw_grid();

    void draw_title();
    void draw_xlabel();
    void draw_ylabel();


    void draw_all_plots();

    PdfDoc *doc;
    PdfPage *page;
    PdfContents *canvas;
    std::string FileName;
  
    // ------------ Paper settings -----------------------
    double PageWidth, PageHeight; // Page size in pts
    double AxisWidth, AxisHeight; // Axis size in pts
    double xTickLength, yTickLength; // tick-length in procentage of axis
    double AxisLineWidth; // Line width of Axis and ticks
    double x_origo, y_origo; // origo position in pts
    double GridLineWidth;
    LineStyle GridLineStyle;

    // ----- Default drawing styles, widths, etc. -------------------
    std::string dFontName;
    int dFontSize;
    double dLineWidth;
    LineStyle dLineStyle;
    
    double MarkerSize; // size of markers in percentage of axis

    // ---- Current Axis settings -----------------
    vec x_ticks, y_ticks; // position of x- and y-ticks in coordinate space
    Array<std::string> xtick_labels, ytick_labels; // label of x- and y-ticks
    double XMin, XMax, YMin, YMax; // axis settings
    double dx, dy; // scale for coordinates to pts
    bool grid; // true if grid should be drawn

    std::string XLabel, YLabel, Title;

    //------ Plot details -------------------------
    int Nplots;
    Array<vec> x_values;
    Array<vec> y_values;
    vec line_width;
    Array<LineStyle> line_style;
    Array<MarkerType> line_marker;
    ivec line_color;
  };

} //namespace itpp

#endif // HAVE_HARU

#endif // __plot_pdf_h
