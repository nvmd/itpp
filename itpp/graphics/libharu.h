/*
 * << H a r u -- Free PDF Library >> -- libharu.h
 *
 * Copyright (c) 1999-2003 Takeshi Kanno <takeshi_kanno@est.hi-ho.ne.jp>
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 * It is provided "as is" without express or implied warranty.
 *
 */

// This file is taken from libhary 1.2.0 beta 3
#define NOPNG
#define NOJPEG
#define NOZLIB


#ifndef _LIB_HARU_H 
#define _LIB_HARU_H 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
#include <exception>

#ifdef _BCC32
#include <exceptio.h>
#endif /* _BCC32 */
#endif /* __cplusplus */

/*----------------------------------------------------------------------------*/
/*------ typedef elements which depending on platform ------------------------*/

typedef int             PdfOID;
typedef char            PdfEntryType;
typedef unsigned char   pdf_uint8;
typedef unsigned short  pdf_uint16;
typedef unsigned long   pdf_uint32;

#ifndef FLT_MIN
#define PDF_PARAM_NODEF 1.17549435e-38F
#else
#define PDF_PARAM_NODEF FLT_MIN
#endif /* FLT_MIN */

#ifndef INT_MIN
#define PDF_INT_MIN     -2147483647
#else
#define PDF_INT_MIN     INT_MIN
#endif /* INT_MIN */

/*----------------------------------------------------------------------------*/
/*---- Using external libraries or not. --------------------------------------*/

#ifdef NOZLIB
#define PDF_ZLIB_MODE           0
#else
#define PDF_ZLIB_MODE           1
#include <zlib.h>
#endif /* NOZLIB */

#ifndef NOPNG
#include <png.h>
#endif /* NOPNG */

#ifndef NOJPEG
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#include <jpeglib.h>
#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* NOJPEG */

/*----------------------------------------------------------------------------*/
/*---- error number definition -----------------------------------------------*/

#define PDF_ERR_MALLOC                  -1
#define PDF_ERR_NEW_OBJECT              -2
#define PDF_ERR_FILE_OPEN               -3
#define PDF_ERR_INVALID_RANGE           -4
#define PDF_RUNTIME_ERROR               -5
#define PDF_ERR_OUT_OF_RANGE            -6
#define PDF_INTERNAL_ERROR              -7
#define PDF_DEFLATEER_ERROR             -8
#define PDF_INVALID_PARAMETER           -9
#define PDF_ERR_LIST_ADD                -10
#define PDF_ERR_INVALID_OBJECT          -11
#define PDF_ERR_INVALID_FONT_WIDTHS     -12
#define PDF_ERR_INVALID_OPERATION       -13
#define PDF_ERR_INVALID_GRAPHIC_MODE    -14
#define PDF_ERR_INVALID_WMF_FILE        -15
#define PDF_ERR_INVALID_CMAP_OBJECT     -16
#define PDF_ERR_NOT_SUPPORTED_FUNCTION  -17
#define PDF_UNKNOWN_ERROR           PDF_INT_MIN

/*----------------------------------------------------------------------------*/
/*------ PreDefined page size ------------------------------------------------*/

#define PDF_PAGE_WIDTH_A4       596
#define PDF_PAGE_HEIGHT_A4      842

#define PDF_DEFAULT_PAGE_WIDTH  PDF_PAGE_WIDTH_A4
#define PDF_DEFAULT_PAGE_HEIGHT PDF_PAGE_HEIGHT_A4

/*----------------------------------------------------------------------------*/
/*------ default values ------------------------------------------------------*/

#define PDF_DEF_BUF_SIZE        1024
#define PDF_DEFLATOR_BUF_SIZE   8192

/*----------------------------------------------------------------------------*/
/*------ collection of flags defining various characteristics of the font. ---*/

#define PDF_FONT_FIXED_WIDTH    1
#define PDF_FONT_SERIF          2
#define PDF_FONT_SYMBOLIC       4
#define PDF_FONT_SCRIPT         8
  /* Reserved                   16 */
#define PDF_FONT_STD_CHARSET    32
#define PDF_FONT_ITALIC         64
  /* Reserved                   128
     Reserved                   256
     Reserved                   512
     Reserved                   1024
     Reserved                   2048
     Reserved                   4096
     Reserved                   8192
     Reserved                   16384
     Reserved                   32768 */
#define PDF_FONT_ALL_CAP        65536
#define PDF_FONT_SMALL_CAP      131072
#define PDF_FONT_FOURCE_BOLD    262144

/*----------------------------------------------------------------------------*/
/*------ limitation of object implementation ---------------------------------*/

#define PDF_LIMIT_MAX_ARRAY         8191
#define PDF_LIMIT_MAX_DICT_ELEMENT  4095
#define PDF_LIMIT_STRING_MAX        65535
#define PDF_LIMIT_MAX_NAME          127
#define PDF_LIMIT_MAX_XREF          250000
#define PDF_LIMIT_MAX_REAL          32767
#define PDF_LIMIT_MIN_REAL          -32767

/*----------------------------------------------------------------------------*/
/*------ limitation of various properties ------------------------------------*/

#define PDF_MIN_HORIZONTALSCALING   10
#define PDF_MAX_HORIZONTALSCALING   300
#define PDF_MIN_WORDSPACE           -30
#define PDF_MAX_WORDSPACE           300
#define PDF_MIN_CHARSPACE           -30
#define PDF_MAX_CHARSPACE           300
#define PDF_MAX_FONTSIZE            300
#define PDF_MAX_ZOOMSIZE            10
#define PDF_MAX_LEADING             300
#define PDF_MAX_LINEWIDTH           100

/*----------------------------------------------------------------------------*/
/*------ default values ------------------------------------------------------*/

#define PDF_DEF_FONT_SIZE           10
#define PDF_DEF_FONT                "Hervetica"
#define PDF_DEF_PAGE_LAYOUT         PDF_SINGLE_PAGE
#define PDF_DEF_PAGE_MODE           PDF_USE_NONE
#define PDF_DEF_WORDSPACE           0
#define PDF_DEF_CHARSPACE           0
#define PDF_DEF_FONTSIZE            10
#define PDF_DEF_HSCALING            100
#define PDF_DEF_LEADING             0
#define PDF_DEF_RENDERING_MODE      PDF_FILL
#define PDF_DEF_RAISE               0
#define PDF_DEF_LINEWIDTH           1
#define PDF_DEF_LINECAP             PDF_BUTT_END
#define PDF_DEF_LINEJOIN            PDF_MITER_JOIN
#define PDF_DEF_MITERLIMIT          10
#define PDF_DEF_PAGE_NUM            1

#define LIB_HARU_VERSION_TXT        "HARU Free PDF Library Version 1.2.0 Beta2"
#define PDF_ERROR_MSG_MAX           512
#define PDF_DEFAULT_ITEMS_PER_BLOCK 20
#define PDF_UNICODE_HEADER_LEN      2

/*----------------------------------------------------------------------------*/
/*------ The pagemode determines. --------------------------------------------*/

enum pdf_page_mode_enum {
    PDF_USE_NONE = 0,
    PDF_USE_OUTLINES,
    PDF_USE_THUMBS,
    PDF_FULL_SCREEN
};
typedef enum pdf_page_mode_enum pdf_page_mode;

/*----------------------------------------------------------------------------*/
/*------ The line cap style --------------------------------------------------*/

enum pdf_line_cap_style_enum {
    PDF_BUTT_END,
    PDF_ROUND_END,
    PDF_PROJECTING_SCUARE_END
};
typedef enum pdf_line_cap_style_enum pdf_line_cap_style;

/*----------------------------------------------------------------------------*/
/*------ The line join style -------------------------------------------------*/

enum pdf_line_join_style_enum {
    PDF_MITER_JOIN = 0,
    PDF_ROUND_JOIN,
    PDF_BEVEL_JOIN
};
typedef enum pdf_line_join_style_enum pdf_line_join_style;

/*----------------------------------------------------------------------------*/
/*------ The text rendering mode ---------------------------------------------*/

enum pdf_text_rendering_mode_enum {
    PDF_FILL = 0,
    PDF_STROKE,
    PDF_FILL_THEN_STROKE,
    PDF_INVISIBLE,
    PDF_FILL_CLIPPING,
    PDF_STROKE_CLIPPING,
    PDF_FILL_STROKE_CLIPPING,
    PDF_CLIPPING,
    PDF_RENDERING_MODE_EOF
};
typedef enum pdf_text_rendering_mode_enum pdf_text_rendering_mode;

/*----------------------------------------------------------------------------*/
/*------ border stype --------------------------------------------------------*/

enum pdf_bs_subtype_enum {
    PDF_BS_SOLID,
    PDF_BS_DASHED,
    PDF_BS_BEVELED,
    PDF_BS_INSET,
    PDF_BS_UNDERLINED
};
typedef enum pdf_bs_subtype_enum pdf_bs_subtype;

/*----------------------------------------------------------------------------*/
/*------ The annotation types ------------------------------------------------*/

enum pdf_annot_type_enum {
    PDF_ANNOT_TEXT_NOTES,
    PDF_ANNOT_LINK,
    PDF_ANNOT_SOUND,
    PDF_ANNOT_FREE_TEXT,
    PDF_ANNOT_STAMP,
    PDF_ANNOT_SQUARE,
    PDF_ANNOT_CIRCLE,
    PDF_ANNOT_STRIKE_OUT,
    PDF_ANNOT_HIGHTLIGHT,
    PDF_ANNOT_UNDERLINE,
    PDF_ANNOT_INK,
    PDF_ANNOT_FILE_ATTACHMENT,
    PDF_ANNOT_POPUP
};
typedef enum pdf_annot_type_enum pdf_annot_type;

enum pdf_annot_flgs_enum {
    PDF_ANNOT_INVISIBLE,
    PDF_ANNOT_HIDDEN,
    PDF_ANNOT_PRINT,
    PDF_ANNOT_NOZOOM,
    PDF_ANNOT_NOROTATE,
    PDF_ANNOT_NOVIEW,
    PDF_ANNOT_READONLY
};
typedef enum pdf_annot_flgs_enum pdf_annot_flgs;

enum pdf_annot_hl_mode_enum {
    PDF_ANNOT_NO_HIGHTLIGHT,
    PDF_ANNOT_INVERT_BOX,
    PDF_ANNOT_INVERT_BORDER,
    PDF_ANNOT_DOWN_APPEARANCE,
    PDF_ANNOT_HL_EOF
};
typedef enum pdf_annot_hl_mode_enum pdf_annot_hl_mode;

enum pdf_annot_icon_names_enum {
    PDF_ANNOT_ICON_COMMENT,
    PDF_ANNOT_ICON_KEY,
    PDF_ANNOT_ICON_NOTE,
    PDF_ANNOT_ICON_HELP,
    PDF_ANNOT_ICON_NEW_PARAGRAPH,
    PDF_ANNOT_ICON_PARAGRAPH,
    PDF_ANNOT_ICON_INSERT,
    PDF_ANNOT_ICON_EOF
};
typedef enum pdf_annot_icon_names_enum pdf_annot_icon_names;

/*----------------------------------------------------------------------------*/
/*------ pdf_destination_type ------------------------------------------------*/

enum pdf_destination_type_enum {
    PDF_XYZ = 0,
    PDF_FIT,
    PDF_FIT_H,
    PDF_FIT_V,
    PDF_FIT_R,
    PDF_FIT_B,
    PDF_FIT_BH,
    PDF_FIT_BV,
    PDF_DST_EOF
};
typedef enum pdf_destination_type_enum pdf_destination_type;

/*----------------------------------------------------------------------------*/
/*------ pdf_page_layout -----------------------------------------------------*/

enum pdf_page_layout_enum {
    PDF_SINGLE_PAGE,
    PDF_ONE_COLUMN,
    PDF_TWO_COLUMN_LEFT,
    PDF_TWO_COLUMN_RIGHT
};
typedef enum pdf_page_layout_enum pdf_page_layout;

/*----------------------------------------------------------------------------*/
/*------ viewer preferences definitions --------------------------------------*/

#define PDF_HIDE_TOOLBAR    1
#define PDF_HIDE_MENUBAR    2
#define PDF_HIDE_WINDOW_UI  4
#define PDF_FIT_WINDOW      8
#define PDF_CENTER_WINDOW   16

/*----------------------------------------------------------------------------*/
/*------ pdf_page_num_style definition ---------------------------------------*/

enum pdf_page_num_style_enum {
    PDF_PAGE_NUM_DECIMAL = 0,
    PDF_PAGE_NUM_UPPER_ROMAN,
    PDF_PAGE_NUM_LOWER_ROMAN,
    PDF_PAGE_NUM_UPPER_LETTERS,
    PDF_PAGE_NUM_LOWER_LETTERS
};
typedef enum pdf_page_num_style_enum pdf_page_num_style;

/*----------------------------------------------------------------------------*/
/*------ pdf_color_space definition ------------------------------------------*/

enum pdf_color_space_enum {
    PDF_CS_DEVICE_GRAY = 0,
    PDF_CS_DEVICE_RGB,
    PDF_CS_DEVICE_CMYK,
    PDF_CS_CAL_GRAY,
    PDF_CS_CAL_RGB,
    PDF_CS_LAB,
    PDF_CS_ICC_BASED,
    PDF_CS_SEPARATION,
    PDF_CS_DEVICE_N,
    PDF_CS_INDEXED,
    PDF_CS_PATTERN,
    PDF_CS_EOF
};
typedef enum pdf_color_space_enum pdf_color_space;

/*----------------------------------------------------------------------------*/
/*------ bits-per-component definition ---------------------------------------*/

#define PDF_BITS_1      1
#define PDF_BITS_2      2
#define PDF_BITS_4      4
#define PDF_BITS_8      8

/*----------------------------------------------------------------------------*/
/*------ predefined font encoding --------------------------------------------*/

enum pdf_encoding_enum {
    PDF_STANDARD_ENCODING = 0,
    PDF_MAC_ROMAN_ENCODING,
    PDF_WIN_ANSI_ENCODING,
    PDF_FONT_SPECIFIC,
    PDF_ENCODING_EOF
};
typedef enum pdf_encoding_enum pdf_encoding;

/*----------------------------------------------------------------------------*/
/*----- Writing Mode ---------------------------------------------------------*/

enum pdf_writing_mode_enum {
    PDF_WMODE_HORIZONTAL = 0,
    PDF_WMODE_VERTICAL
};
typedef enum pdf_writing_mode_enum pdf_writing_mode;

/*----------------------------------------------------------------------------*/
/*----- Graphics Mode --------------------------------------------------------*/

enum pdf_graphics_mode_enum {
    PDF_GMODE_PAGE_DESCRIPTION = 0,
    PDF_GMODE_PATH_OBJECT,
    PDF_GMODE_TEXT_OBJECT,
    PDF_GMODE_CLIPPING_PATH,
    PDF_GMODE_SHADING,
    PDF_GMODE_INLINE_IMAGE,
    PDF_GMODE_EXTERNAL_OBJECT,
    PDF_GMODE_EOF
};
typedef enum pdf_graphics_mode_enum pdf_graphics_mode;

/*----------------------------------------------------------------------------*/
/*----- filter type definition -----------------------------------------------*/

#define PDF_FILTER_NONE         0
#define PDF_FILTER_DEFLATE      2
#define PDF_FILTER_DCT_DECODE   8

typedef unsigned int      pdf_filter;

/*----------------------------------------------------------------------------*/
/*----- proc set element definition ------------------------------------------*/

#define PDF_PROCSET_PDF         2
#define PDF_PROCSET_TEXT        4
#define PDF_PROCSET_IMAGEB      8
#define PDF_PROCSET_IMAGEC      16
#define PDF_PROCSET_IMAGEI      32

/*----------------------------------------------------------------------------*/
/*----- pdf_point struct -----------------------------------------------------*/

struct pdf_point_struct {
        double x;
        double y;
};
typedef struct pdf_point_struct pdf_point;

/*----------------------------------------------------------------------------*/
/*----- pdf_rect struct ------------------------------------------------------*/

struct pdf_rect_struct {
        double left;
        double bottom;
        double right;
        double top;
};
typedef struct pdf_rect_struct pdf_rect;

/*----------------------------------------------------------------------------*/
/*----- pdf_box struct -------------------------------------------------------*/

struct pdf_box_struct {
        int left;
        int bottom;
        int right;
        int top;
};
typedef struct pdf_box_struct pdf_box;

/*----------------------------------------------------------------------------*/
/*----- pdf_text_matrix struct -----------------------------------------------*/

struct pdf_text_matrix_struct {
        double a;
        double b;
        double c;
        double d;
        double x;
        double y;
};
typedef struct pdf_text_matrix_struct pdf_text_matrix;

/*----------------------------------------------------------------------------*/
/*----- pdf_text_width struct ------------------------------------------------*/

struct pdf_text_width_struct {
        int numchars;
        int numwords;
        int width;
};
typedef struct pdf_text_width_struct pdf_text_width;

/*----------------------------------------------------------------------------*/
/*----- pdf_date struct ------------------------------------------------------*/

struct pdf_date_struct {
        int year;
        int month;
        int day;
        int hour;
        int minutes;
        int seconds;
        char ind;
        int off_hour;
        int off_minutes;
};
typedef struct pdf_date_struct pdf_date;

/*----------------------------------------------------------------------------*/
/*----- pdf_rgb_color struct -------------------------------------------------*/

struct pdf_pal_color_struct {
        unsigned char red;
        unsigned char green;
        unsigned char blue;
};
typedef struct pdf_pal_color_struct pdf_pal_color;

struct pdf_rgb_color_struct {
        double red;
        double green;
        double blue;
};
typedef struct pdf_rgb_color_struct pdf_rgb_color;

/*---------------------------------------------------------------------------*/
/*----- definition for multibyte fonts --------------------------------------*/

typedef unsigned short pdf_cid;

typedef struct pdf_cid_range_struct {
    unsigned short from;
    unsigned short to;
    pdf_cid cid;
} pdf_cid_range;

typedef struct pdf_mb_unicode_struct1 {
    unsigned short mbchar;
    unsigned short unicode;
} pdf_mb_unicode_map1;

typedef struct pdf_mb_unicode_struct2 {
    unsigned short from;
    unsigned short to;
    unsigned short unicode;
} pdf_mb_unicode_map2;

typedef enum pdf_cid_widths_type_enum {
    PDF_CID_W_TYPE_FROM_TO = 0,
    PDF_CID_W_TYPE_FROM_ARRAY
} pdf_cid_widths_type;

typedef struct pdf_cid_width_struct {
    pdf_cid_widths_type type;
    pdf_cid from_cid;
    pdf_cid to_cid;
    unsigned int* widths;
} pdf_cid_width;

typedef enum pdf_byte_type_enum {
    PDF_BYTE_TYPE_SINGLE = 0,
    PDF_BYTE_TYPE_LEAD,
    PDF_BYTE_TYPE_TRIAL
} pdf_byte_type;

/*----------------------------------------------------------------------------*/
/*----- definition for font encoding -----------------------------------------*/

#define PDF_CHAR_NOTDEF ".notdef" 

/*----------------------------------------------------------------------------*/
/*----- country code definition ----------------------------------------------*/

#define PDF_COUNTRY_AF  "AF"    /* AFGHANISTAN */
#define PDF_COUNTRY_AL  "AL"    /* ALBANIA */
#define PDF_COUNTRY_DZ  "DZ"    /* ALGERIA */
#define PDF_COUNTRY_AS  "AS"    /* AMERICAN SAMOA */
#define PDF_COUNTRY_AD  "AD"    /* ANDORRA */
#define PDF_COUNTRY_AO  "AO"    /* ANGOLA */
#define PDF_COUNTRY_AI  "AI"    /* ANGUILLA */
#define PDF_COUNTRY_AQ  "AQ"    /* ANTARCTICA */
#define PDF_COUNTRY_AG  "AG"    /* ANTIGUA AND BARBUDA */
#define PDF_COUNTRY_AR  "AR"    /* ARGENTINA */
#define PDF_COUNTRY_AM  "AM"    /* ARMENIA */
#define PDF_COUNTRY_AW  "AW"    /* ARUBA */
#define PDF_COUNTRY_AU  "AU"    /* AUSTRALIA */
#define PDF_COUNTRY_AT  "AT"    /* AUSTRIA */
#define PDF_COUNTRY_AZ  "AZ"    /* AZERBAIJAN */
#define PDF_COUNTRY_BS  "BS"    /* BAHAMAS */
#define PDF_COUNTRY_BH  "BH"    /* BAHRAIN */
#define PDF_COUNTRY_BD  "BD"    /* BANGLADESH */
#define PDF_COUNTRY_BB  "BB"    /* BARBADOS */
#define PDF_COUNTRY_BY  "BY"    /* BELARUS */
#define PDF_COUNTRY_BE  "BE"    /* BELGIUM */
#define PDF_COUNTRY_BZ  "BZ"    /* BELIZE */
#define PDF_COUNTRY_BJ  "BJ"    /* BENIN */
#define PDF_COUNTRY_BM  "BM"    /* BERMUDA */
#define PDF_COUNTRY_BT  "BT"    /* BHUTAN */
#define PDF_COUNTRY_BO  "BO"    /* BOLIVIA */
#define PDF_COUNTRY_BA  "BA"    /* BOSNIA AND HERZEGOWINA */
#define PDF_COUNTRY_BW  "BW"    /* BOTSWANA */
#define PDF_COUNTRY_BV  "BV"    /* BOUVET ISLAND */
#define PDF_COUNTRY_BR  "BR"    /* BRAZIL */
#define PDF_COUNTRY_IO  "IO"    /* BRITISH INDIAN OCEAN TERRITORY */
#define PDF_COUNTRY_BN  "BN"    /* BRUNEI DARUSSALAM */
#define PDF_COUNTRY_BG  "BG"    /* BULGARIA */
#define PDF_COUNTRY_BF  "BF"    /* BURKINA FASO */
#define PDF_COUNTRY_BI  "BI"    /* BURUNDI */
#define PDF_COUNTRY_KH  "KH"    /* CAMBODIA */
#define PDF_COUNTRY_CM  "CM"    /* CAMEROON */
#define PDF_COUNTRY_CA  "CA"    /* CANADA */
#define PDF_COUNTRY_CV  "CV"    /* CAPE VERDE */
#define PDF_COUNTRY_KY  "KY"    /* CAYMAN ISLANDS */
#define PDF_COUNTRY_CF  "CF"    /* CENTRAL AFRICAN REPUBLIC */
#define PDF_COUNTRY_TD  "TD"    /* CHAD */
#define PDF_COUNTRY_CL  "CL"    /* CHILE */
#define PDF_COUNTRY_CN  "CN"    /* CHINA */
#define PDF_COUNTRY_CX  "CX"    /* CHRISTMAS ISLAND */
#define PDF_COUNTRY_CC  "CC"    /* COCOS (KEELING) ISLANDS */
#define PDF_COUNTRY_CO  "CO"    /* COLOMBIA */
#define PDF_COUNTRY_KM  "KM"    /* COMOROS */
#define PDF_COUNTRY_CG  "CG"    /* CONGO */
#define PDF_COUNTRY_CK  "CK"    /* COOK ISLANDS */
#define PDF_COUNTRY_CR  "CR"    /* COSTA RICA */
#define PDF_COUNTRY_CI  "CI"    /* COTE D'IVOIRE */
#define PDF_COUNTRY_HR  "HR"    /* CROATIA (local name: Hrvatska) */
#define PDF_COUNTRY_CU  "CU"    /* CUBA */
#define PDF_COUNTRY_CY  "CY"    /* CYPRUS */
#define PDF_COUNTRY_CZ  "CZ"    /* CZECH REPUBLIC */
#define PDF_COUNTRY_DK  "DK"    /* DENMARK */
#define PDF_COUNTRY_DJ  "DJ"    /* DJIBOUTI */
#define PDF_COUNTRY_DM  "DM"    /* DOMINICA */
#define PDF_COUNTRY_DO  "DO"    /* DOMINICAN REPUBLIC */
#define PDF_COUNTRY_TP  "TP"    /* EAST TIMOR */
#define PDF_COUNTRY_EC  "EC"    /* ECUADOR */
#define PDF_COUNTRY_EG  "EG"    /* EGYPT */
#define PDF_COUNTRY_SV  "SV"    /* EL SALVADOR */
#define PDF_COUNTRY_GQ  "GQ"    /* EQUATORIAL GUINEA */
#define PDF_COUNTRY_ER  "ER"    /* ERITREA */
#define PDF_COUNTRY_EE  "EE"    /* ESTONIA */
#define PDF_COUNTRY_ET  "ET"    /* ETHIOPIA */
#define PDF_COUNTRY_FK  "FK"   /* FALKLAND ISLANDS (MALVINAS) */
#define PDF_COUNTRY_FO  "FO"    /* FAROE ISLANDS */
#define PDF_COUNTRY_FJ  "FJ"    /* FIJI */
#define PDF_COUNTRY_FI  "FI"    /* FINLAND */
#define PDF_COUNTRY_FR  "FR"    /* FRANCE */
#define PDF_COUNTRY_FX  "FX"    /* FRANCE, METROPOLITAN */
#define PDF_COUNTRY_GF  "GF"    /* FRENCH GUIANA */
#define PDF_COUNTRY_PF  "PF"    /* FRENCH POLYNESIA */
#define PDF_COUNTRY_TF  "TF"    /* FRENCH SOUTHERN TERRITORIES */
#define PDF_COUNTRY_GA  "GA"    /* GABON */
#define PDF_COUNTRY_GM  "GM"    /* GAMBIA */
#define PDF_COUNTRY_GE  "GE"    /* GEORGIA */
#define PDF_COUNTRY_DE  "DE"    /* GERMANY */
#define PDF_COUNTRY_GH  "GH"    /* GHANA */
#define PDF_COUNTRY_GI  "GI"    /* GIBRALTAR */
#define PDF_COUNTRY_GR  "GR"    /* GREECE */
#define PDF_COUNTRY_GL  "GL"    /* GREENLAND */
#define PDF_COUNTRY_GD  "GD"    /* GRENADA */
#define PDF_COUNTRY_GP  "GP"    /* GUADELOUPE */
#define PDF_COUNTRY_GU  "GU"    /* GUAM */
#define PDF_COUNTRY_GT  "GT"    /* GUATEMALA */
#define PDF_COUNTRY_GN  "GN"    /* GUINEA */
#define PDF_COUNTRY_GW  "GW"    /* GUINEA-BISSAU */
#define PDF_COUNTRY_GY  "GY"    /* GUYANA */
#define PDF_COUNTRY_HT  "HT"    /* HAITI */
#define PDF_COUNTRY_HM  "HM"    /* HEARD AND MC DONALD ISLANDS */
#define PDF_COUNTRY_HN  "HN"    /* HONDURAS */
#define PDF_COUNTRY_HK  "HK"    /* HONG KONG */
#define PDF_COUNTRY_HU  "HU"    /* HUNGARY */
#define PDF_COUNTRY_IS  "IS"    /* ICELAND */
#define PDF_COUNTRY_IN  "IN"    /* INDIA */
#define PDF_COUNTRY_ID  "ID"    /* INDONESIA */
#define PDF_COUNTRY_IR  "IR"    /* IRAN (ISLAMIC REPUBLIC OF) */
#define PDF_COUNTRY_IQ  "IQ"    /* IRAQ */
#define PDF_COUNTRY_IE  "IE"    /* IRELAND */
#define PDF_COUNTRY_IL  "IL"    /* ISRAEL */
#define PDF_COUNTRY_IT  "IT"    /* ITALY */
#define PDF_COUNTRY_JM  "JM"    /* JAMAICA */
#define PDF_COUNTRY_JP  "JP"    /* JAPAN */
#define PDF_COUNTRY_JO  "JO"    /* JORDAN */
#define PDF_COUNTRY_KZ  "KZ"    /* KAZAKHSTAN */
#define PDF_COUNTRY_KE  "KE"    /* KENYA */
#define PDF_COUNTRY_KI  "KI"    /* KIRIBATI */
#define PDF_COUNTRY_KP  "KP"    /* KOREA, DEMOCRATIC PEOPLE'S REPUBLIC OF */
#define PDF_COUNTRY_KR  "KR"    /* KOREA, REPUBLIC OF */
#define PDF_COUNTRY_KW  "KW"    /* KUWAIT */
#define PDF_COUNTRY_KG  "KG"    /* KYRGYZSTAN */
#define PDF_COUNTRY_LA  "LA"    /* LAO PEOPLE'S DEMOCRATIC REPUBLIC */
#define PDF_COUNTRY_LV  "LV"    /* LATVIA */
#define PDF_COUNTRY_LB  "LB"    /* LEBANON */
#define PDF_COUNTRY_LS  "LS"    /* LESOTHO */
#define PDF_COUNTRY_LR  "LR"    /* LIBERIA */
#define PDF_COUNTRY_LY  "LY"    /* LIBYAN ARAB JAMAHIRIYA */
#define PDF_COUNTRY_LI  "LI"    /* LIECHTENSTEIN */
#define PDF_COUNTRY_LT  "LT"    /* LITHUANIA */
#define PDF_COUNTRY_LU  "LU"    /* LUXEMBOURG */
#define PDF_COUNTRY_MO  "MO"    /* MACAU */
#define PDF_COUNTRY_MK  "MK"   /* MACEDONIA, THE FORMER YUGOSLAV REPUBLIC OF */
#define PDF_COUNTRY_MG  "MG"    /* MADAGASCAR */
#define PDF_COUNTRY_MW  "MW"    /* MALAWI */
#define PDF_COUNTRY_MY  "MY"    /* MALAYSIA */
#define PDF_COUNTRY_MV  "MV"    /* MALDIVES */
#define PDF_COUNTRY_ML  "ML"    /* MALI */
#define PDF_COUNTRY_MT  "MT"    /* MALTA */
#define PDF_COUNTRY_MH  "MH"    /* MARSHALL ISLANDS */
#define PDF_COUNTRY_MQ  "MQ"    /* MARTINIQUE */
#define PDF_COUNTRY_MR  "MR"    /* MAURITANIA */
#define PDF_COUNTRY_MU  "MU"    /* MAURITIUS */
#define PDF_COUNTRY_YT  "YT"    /* MAYOTTE */
#define PDF_COUNTRY_MX  "MX"    /* MEXICO */
#define PDF_COUNTRY_FM  "FM"    /* MICRONESIA, FEDERATED STATES OF */
#define PDF_COUNTRY_MD  "MD"    /* MOLDOVA, REPUBLIC OF */
#define PDF_COUNTRY_MC  "MC"    /* MONACO */
#define PDF_COUNTRY_MN  "MN"    /* MONGOLIA */
#define PDF_COUNTRY_MS  "MS"    /* MONTSERRAT */
#define PDF_COUNTRY_MA  "MA"    /* MOROCCO */
#define PDF_COUNTRY_MZ  "MZ"    /* MOZAMBIQUE */
#define PDF_COUNTRY_MM  "MM"    /* MYANMAR */
#define PDF_COUNTRY_NA  "NA"    /* NAMIBIA */
#define PDF_COUNTRY_NR  "NR"    /* NAURU */
#define PDF_COUNTRY_NP  "NP"    /* NEPAL */
#define PDF_COUNTRY_NL  "NL"    /* NETHERLANDS */
#define PDF_COUNTRY_AN  "AN"    /* NETHERLANDS ANTILLES */
#define PDF_COUNTRY_NC  "NC"    /* NEW CALEDONIA */
#define PDF_COUNTRY_NZ  "NZ"    /* NEW ZEALAND */
#define PDF_COUNTRY_NI  "NI"    /* NICARAGUA */
#define PDF_COUNTRY_NE  "NE"    /* NIGER */
#define PDF_COUNTRY_NG  "NG"    /* NIGERIA */
#define PDF_COUNTRY_NU  "NU"    /* NIUE */
#define PDF_COUNTRY_NF  "NF"    /* NORFOLK ISLAND */
#define PDF_COUNTRY_MP  "MP"    /* NORTHERN MARIANA ISLANDS */
#define PDF_COUNTRY_NO  "NO"    /* NORWAY */
#define PDF_COUNTRY_OM  "OM"    /* OMAN */
#define PDF_COUNTRY_PK  "PK"    /* PAKISTAN */
#define PDF_COUNTRY_PW  "PW"    /* PALAU */
#define PDF_COUNTRY_PA  "PA"    /* PANAMA */
#define PDF_COUNTRY_PG  "PG"    /* PAPUA NEW GUINEA */
#define PDF_COUNTRY_PY  "PY"    /* PARAGUAY */
#define PDF_COUNTRY_PE  "PE"    /* PERU */
#define PDF_COUNTRY_PH  "PH"    /* PHILIPPINES */
#define PDF_COUNTRY_PN  "PN"    /* PITCAIRN */
#define PDF_COUNTRY_PL  "PL"    /* POLAND */
#define PDF_COUNTRY_PT  "PT"    /* PORTUGAL */
#define PDF_COUNTRY_PR  "PR"    /* PUERTO RICO */
#define PDF_COUNTRY_QA  "QA"    /* QATAR */
#define PDF_COUNTRY_RE  "RE"    /* REUNION */
#define PDF_COUNTRY_RO  "RO"    /* ROMANIA */
#define PDF_COUNTRY_RU  "RU"    /* RUSSIAN FEDERATION */
#define PDF_COUNTRY_RW  "RW"    /* RWANDA */
#define PDF_COUNTRY_KN  "KN"    /* SAINT KITTS AND NEVIS */
#define PDF_COUNTRY_LC  "LC"    /* SAINT LUCIA */
#define PDF_COUNTRY_VC  "VC"    /* SAINT VINCENT AND THE GRENADINES */
#define PDF_COUNTRY_WS  "WS"    /* SAMOA */
#define PDF_COUNTRY_SM  "SM"    /* SAN MARINO */
#define PDF_COUNTRY_ST  "ST"    /* SAO TOME AND PRINCIPE */
#define PDF_COUNTRY_SA  "SA"    /* SAUDI ARABIA */
#define PDF_COUNTRY_SN  "SN"    /* SENEGAL */
#define PDF_COUNTRY_SC  "SC"    /* SEYCHELLES */
#define PDF_COUNTRY_SL  "SL"    /* SIERRA LEONE */
#define PDF_COUNTRY_SG  "SG"    /* SINGAPORE */
#define PDF_COUNTRY_SK  "SK"    /* SLOVAKIA (Slovak Republic) */
#define PDF_COUNTRY_SI  "SI"    /* SLOVENIA */
#define PDF_COUNTRY_SB  "SB"    /* SOLOMON ISLANDS */
#define PDF_COUNTRY_SO  "SO"    /* SOMALIA */
#define PDF_COUNTRY_ZA  "ZA"    /* SOUTH AFRICA */
#define PDF_COUNTRY_ES  "ES"    /* SPAIN */
#define PDF_COUNTRY_LK  "LK"    /* SRI LANKA */
#define PDF_COUNTRY_SH  "SH"    /* ST. HELENA */
#define PDF_COUNTRY_PM  "PM"    /* ST. PIERRE AND MIQUELON */
#define PDF_COUNTRY_SD  "SD"    /* SUDAN */
#define PDF_COUNTRY_SR  "SR"    /* SURINAME */
#define PDF_COUNTRY_SJ  "SJ"    /* SVALBARD AND JAN MAYEN ISLANDS */
#define PDF_COUNTRY_SZ  "SZ"    /* SWAZILAND */
#define PDF_COUNTRY_SE  "SE"    /* SWEDEN */
#define PDF_COUNTRY_CH  "CH"    /* SWITZERLAND */
#define PDF_COUNTRY_SY  "SY"    /* SYRIAN ARAB REPUBLIC */
#define PDF_COUNTRY_TW  "TW"    /* TAIWAN, PROVINCE OF CHINA */
#define PDF_COUNTRY_TJ  "TJ"    /* TAJIKISTAN */
#define PDF_COUNTRY_TZ  "TZ"    /* TANZANIA, UNITED REPUBLIC OF */
#define PDF_COUNTRY_TH  "TH"    /* THAILAND */
#define PDF_COUNTRY_TG  "TG"    /* TOGO */
#define PDF_COUNTRY_TK  "TK"    /* TOKELAU */
#define PDF_COUNTRY_TO  "TO"    /* TONGA */
#define PDF_COUNTRY_TT  "TT"    /* TRINIDAD AND TOBAGO */
#define PDF_COUNTRY_TN  "TN"    /* TUNISIA */
#define PDF_COUNTRY_TR  "TR"    /* TURKEY */
#define PDF_COUNTRY_TM  "TM"    /* TURKMENISTAN */
#define PDF_COUNTRY_TC  "TC"    /* TURKS AND CAICOS ISLANDS */
#define PDF_COUNTRY_TV  "TV"    /* TUVALU */
#define PDF_COUNTRY_UG  "UG"    /* UGANDA */
#define PDF_COUNTRY_UA  "UA"    /* UKRAINE */
#define PDF_COUNTRY_AE  "AE"    /* UNITED ARAB EMIRATES */
#define PDF_COUNTRY_GB  "GB"    /* UNITED KINGDOM */
#define PDF_COUNTRY_US  "US"    /* UNITED STATES */
#define PDF_COUNTRY_UM  "UM"    /* UNITED STATES MINOR OUTLYING ISLANDS */
#define PDF_COUNTRY_UY  "UY"    /* URUGUAY */
#define PDF_COUNTRY_UZ  "UZ"    /* UZBEKISTAN */
#define PDF_COUNTRY_VU  "VU"    /* VANUATU */
#define PDF_COUNTRY_VA  "VA"    /* VATICAN CITY STATE (HOLY SEE) */
#define PDF_COUNTRY_VE  "VE"    /* VENEZUELA */
#define PDF_COUNTRY_VN  "VN"    /* VIET NAM */
#define PDF_COUNTRY_VG  "VG"    /* VIRGIN ISLANDS (BRITISH) */
#define PDF_COUNTRY_VI  "VI"    /* VIRGIN ISLANDS (U.S.) */
#define PDF_COUNTRY_WF  "WF"    /* WALLIS AND FUTUNA ISLANDS */
#define PDF_COUNTRY_EH  "EH"    /* WESTERN SAHARA */
#define PDF_COUNTRY_YE  "YE"    /* YEMEN */
#define PDF_COUNTRY_YU  "YU"    /* YUGOSLAVIA */
#define PDF_COUNTRY_ZR  "ZR"    /* ZAIRE */
#define PDF_COUNTRY_ZM  "ZM"    /* ZAMBIA */
#define PDF_COUNTRY_ZW  "ZW"    /* ZIMBABWE */

/*----------------------------------------------------------------------------*/
/*----- lang code definition -------------------------------------------------*/

#define PDF_LANG_AA    "aa"     /* Afar */
#define PDF_LANG_AB    "ab"     /* Abkhazian */
#define PDF_LANG_AF    "af"     /* Afrikaans */
#define PDF_LANG_AM    "am"     /* Amharic */
#define PDF_LANG_AR    "ar"     /* Arabic */
#define PDF_LANG_AS    "as"     /* Assamese */
#define PDF_LANG_AY    "ay"     /* Aymara */
#define PDF_LANG_AZ    "az"     /* Azerbaijani */
#define PDF_LANG_BA    "ba"     /* Bashkir */
#define PDF_LANG_BE    "be"     /* Byelorussian */
#define PDF_LANG_BG    "bg"     /* Bulgarian */
#define PDF_LANG_BH    "bh"     /* Bihari */
#define PDF_LANG_BI    "bi"     /* Bislama */
#define PDF_LANG_BN    "bn"     /* Bengali Bangla */
#define PDF_LANG_BO    "bo"     /* Tibetan */
#define PDF_LANG_BR    "br"     /* Breton */
#define PDF_LANG_CA    "ca"     /* Catalan */
#define PDF_LANG_CO    "co"     /* Corsican */
#define PDF_LANG_CS    "cs"     /* Czech */
#define PDF_LANG_CY    "cy"     /* Welsh */
#define PDF_LANG_DA    "da"     /* Danish */
#define PDF_LANG_DE    "de"     /* German */
#define PDF_LANG_DZ    "dz"     /* Bhutani */
#define PDF_LANG_EL    "el"     /* Greek */
#define PDF_LANG_EN    "en"     /* English */
#define PDF_LANG_EO    "eo"     /* Esperanto */
#define PDF_LANG_ES    "es"     /* Spanish */
#define PDF_LANG_ET    "et"     /* Estonian */
#define PDF_LANG_EU    "eu"     /* Basque */
#define PDF_LANG_FA    "fa"     /* Persian */
#define PDF_LANG_FI    "fi"     /* Finnish */
#define PDF_LANG_FJ    "fj"     /* Fiji */
#define PDF_LANG_FO    "fo"     /* Faeroese */
#define PDF_LANG_FR    "fr"     /* French */
#define PDF_LANG_FY    "fy"     /* Frisian */
#define PDF_LANG_GA    "ga"     /* Irish */
#define PDF_LANG_GD    "gd"     /* Scots Gaelic */
#define PDF_LANG_GL    "gl"     /* Galician */
#define PDF_LANG_GN    "gn"     /* Guarani */
#define PDF_LANG_GU    "gu"     /* Gujarati */
#define PDF_LANG_HA    "ha"     /* Hausa */
#define PDF_LANG_HI    "hi"     /* Hindi */
#define PDF_LANG_HR    "hr"     /* Croatian */
#define PDF_LANG_HU    "hu"     /* Hungarian */
#define PDF_LANG_HY    "hy"     /* Armenian */
#define PDF_LANG_IA    "ia"     /* Interlingua */
#define PDF_LANG_IE    "ie"     /* Interlingue */
#define PDF_LANG_IK    "ik"     /* Inupiak */
#define PDF_LANG_IN    "in"     /* Indonesian */
#define PDF_LANG_IS    "is"     /* Icelandic */
#define PDF_LANG_IT    "it"     /* Italian */
#define PDF_LANG_IW    "iw"     /* Hebrew */
#define PDF_LANG_JA    "ja"     /* Japanese */
#define PDF_LANG_JI    "ji"     /* Yiddish */
#define PDF_LANG_JW    "jw"     /* Javanese */
#define PDF_LANG_KA    "ka"     /* Georgian */
#define PDF_LANG_KK    "kk"     /* Kazakh */
#define PDF_LANG_KL    "kl"     /* Greenlandic */
#define PDF_LANG_KM    "km"     /* Cambodian */
#define PDF_LANG_KN    "kn"     /* Kannada */
#define PDF_LANG_KO    "ko"     /* Korean */
#define PDF_LANG_KS    "ks"     /* Kashmiri */
#define PDF_LANG_KU    "ku"     /* Kurdish */
#define PDF_LANG_KY    "ky"     /* Kirghiz */
#define PDF_LANG_LA    "la"     /* Latin */
#define PDF_LANG_LN    "ln"     /* Lingala */
#define PDF_LANG_LO    "lo"     /* Laothian */
#define PDF_LANG_LT    "lt"     /* Lithuanian */
#define PDF_LANG_LV    "lv"     /* Latvian,Lettish */
#define PDF_LANG_MG    "mg"     /* Malagasy */
#define PDF_LANG_MI    "mi"     /* Maori */
#define PDF_LANG_MK    "mk"     /* Macedonian */
#define PDF_LANG_ML    "ml"     /* Malayalam */
#define PDF_LANG_MN    "mn"     /* Mongolian */
#define PDF_LANG_MO    "mo"     /* Moldavian */
#define PDF_LANG_MR    "mr"     /* Marathi */
#define PDF_LANG_MS    "ms"     /* Malay */
#define PDF_LANG_MT    "mt"     /* Maltese */
#define PDF_LANG_MY    "my"     /* Burmese */
#define PDF_LANG_NA    "na"     /* Nauru */
#define PDF_LANG_NE    "ne"     /* Nepali */
#define PDF_LANG_NL    "nl"     /* Dutch */
#define PDF_LANG_NO    "no"     /* Norwegian */
#define PDF_LANG_OC    "oc"     /* Occitan */
#define PDF_LANG_OM    "om"     /* (Afan)Oromo */
#define PDF_LANG_OR    "or"     /* Oriya */
#define PDF_LANG_PA    "pa"     /* Punjabi */
#define PDF_LANG_PL    "pl"     /* Polish */
#define PDF_LANG_PS    "ps"     /* Pashto,Pushto */
#define PDF_LANG_PT    "pt"     /* Portuguese  */
#define PDF_LANG_QU    "qu"     /* Quechua */
#define PDF_LANG_RM    "rm"     /* Rhaeto-Romance */
#define PDF_LANG_RN    "rn"     /* Kirundi */
#define PDF_LANG_RO    "ro"     /* Romanian */
#define PDF_LANG_RU    "ru"     /* Russian */
#define PDF_LANG_RW    "rw"     /* Kinyarwanda */
#define PDF_LANG_SA    "sa"     /* Sanskrit */
#define PDF_LANG_SD    "sd"     /* Sindhi */
#define PDF_LANG_SG    "sg"     /* Sangro */
#define PDF_LANG_SH    "sh"     /* Serbo-Croatian */
#define PDF_LANG_SI    "si"     /* Singhalese */
#define PDF_LANG_SK    "sk"     /* Slovak */
#define PDF_LANG_SL    "sl"     /* Slovenian */
#define PDF_LANG_SM    "sm"     /* Samoan */
#define PDF_LANG_SN    "sn"     /* Shona */
#define PDF_LANG_SO    "so"     /* Somali */
#define PDF_LANG_SQ    "sq"     /* Albanian */
#define PDF_LANG_SR    "sr"     /* Serbian */
#define PDF_LANG_SS    "ss"     /* Siswati */
#define PDF_LANG_ST    "st"     /* Sesotho */
#define PDF_LANG_SU    "su"     /* Sundanese */
#define PDF_LANG_SV    "sv"     /* Swedish */
#define PDF_LANG_SW    "sw"     /* Swahili */
#define PDF_LANG_TA    "ta"     /* Tamil */
#define PDF_LANG_TE    "te"     /* Tegulu */
#define PDF_LANG_TG    "tg"     /* Tajik */
#define PDF_LANG_TH    "th"     /* Thai */
#define PDF_LANG_TI    "ti"     /* Tigrinya */
#define PDF_LANG_TK    "tk"     /* Turkmen */
#define PDF_LANG_TL    "tl"     /* Tagalog */
#define PDF_LANG_TN    "tn"     /* Setswanato Tonga */
#define PDF_LANG_TR    "tr"     /* Turkish */
#define PDF_LANG_TS    "ts"     /* Tsonga */
#define PDF_LANG_TT    "tt"     /* Tatar */
#define PDF_LANG_TW    "tw"     /* Twi */
#define PDF_LANG_UK    "uk"     /* Ukrainian */
#define PDF_LANG_UR    "ur"     /* Urdu */
#define PDF_LANG_UZ    "uz"     /* Uzbek */
#define PDF_LANG_VI    "vi"     /* Vietnamese */
#define PDF_LANG_VO    "vo"     /* Volapuk */
#define PDF_LANG_WO    "wo"     /* Wolof */
#define PDF_LANG_XH    "xh"     /* Xhosa */
#define PDF_LANG_YO    "yo"     /* Yoruba */
#define PDF_LANG_ZH    "zh"     /* Chinese */
#define PDF_LANG_ZU    "zu"     /* Zulu */

/*----------------------------------------------------------------------------*/
/*----- permission flags (only Revision 2 is supported)-----------------------*/

#define PDF_ENABLE_READ         0
#define PDF_ENABLE_PRINT        4 
#define PDF_ENABLE_EDIT_ALL     8
#define PDF_ENABLE_COPY         16
#define PDF_ENABLE_EDIT         32

/*----------------------------------------------------------------------------*/
/*----- default values for encrypt -------------------------------------------*/

#define PDF_ID_LEN              16
#define PDF_PASSWD_LEN          32
#define PDF_ENCRYPT_KEY_LEN     5
#define PDF_MD5_KEY_LEN         16
#define PDF_PERMISSION_PAD      0xFFFFFFC0
                             
/*----------------------------------------------------------------------------*/
/*----- check macros ---------------------------------------------------------*/

#define PDF_CHECK_RECT(RECT) \
    (RECT.left < RECT.right && RECT.bottom < RECT.top)

/*----------------------------------------------------------------------------*/
/*----- DEBUG macros ---------------------------------------------------------*/

#ifdef DEBUG
#define PDF_DEBUG_PRINT(ARGS) printf ARGS
#else /* DEBUG == 0 */
#define PDF_DEBUG_PRINT(ARGS) void(0)
#endif

void pdf_print_binary(const unsigned char* buf, int len, char* caption);

#ifdef DEBUG
#define PDF_PRINT_BINARY(BUF, LEN, CAPTION) pdf_print_binary(BUF, LEN, CAPTION)
#else
#define PDF_PRINT_BINARY(BUF, LEN, CAPTION) void(0)
#endif
 
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus

#ifdef _NOT_SUPPORT_STD

#define PDF_STD_EXCEPTION   exception
#define ALLOC_ERROR     exception

int throw_new_handler(size_t size);
void throw_new_handler();

#else /* _NOT_SUPPORT_STD */

#define ALLOC_ERROR         std::bad_alloc
#define PDF_STD_EXCEPTION   std::exception

#endif /* _NOT_SUPPORT_STD */

/*----------------------------------------------------------------------------*/
/*----- auto-ptr object-id ---------------------------------------------------*/

typedef enum pdf_auto_ptr_object_type_enum {
    PDF_OBJECT_UNKNOWN = 0,
    PDF_OBJECT_FONT_DEF,
    PDF_OBJECT_ENCODING_DEF,
    PDF_OBJECT_CMAP
} pdf_auto_ptr_object_type;

/*----------------------------------------------------------------------------*/
/*----- PdfException class ---------------------------------------------------*/

class PdfException : public PDF_STD_EXCEPTION
{
public:
                    PdfException(int code, const char* fmt, ...);
#ifdef _NO_EXCEPT_LIST
                    ~PdfException() {};
        const char* what() const;
#else
                    ~PdfException() throw() {};
        const char* what() const throw();
#endif
        int         GetCode()       { return fCode; }
private:
        char        fErrorBuf[512];
        int         fCode;
};

/*----------------------------------------------------------------------------*/
/*----- PdfList class --------------------------------------------------------*/

class PdfList {

public:
                    PdfList(int itemsPerBlock = PDF_DEFAULT_ITEMS_PER_BLOCK);
                    ~PdfList();

        bool        AddItem(void* item);
        bool        AddItem(void* item, int atIndex);
        bool        RemoveItem(const void* item)
                            { return (RemoveItem(IndexOf(item)) != NULL); }
        void*       RemoveItem(int index);
        void        Clear()             { fCount = 0; }

        void*       ItemAt(int index)   
                       { return (fCount <= index || index < 0) ? NULL : 
                           fObjectList[index]; }
        bool        HasItem(const void* item)
                                        { return (IndexOf(item) > 0); }
        int         IndexOf(const void* item);
        int         CountItems()        { return fCount; }
        bool        IsEmpty()           { return (fCount != 0); }

private:
        int         fBlockSize;
        int         fItemsPerBlock;
        int         fCount;
        void**      fObjectList;

        bool        Resize(int count);
};

/*----------------------------------------------------------------------------*/
/*----- PdfStreamBase --------------------------------------------------------*/

class PdfStreamBase
{
public:
                            PdfStreamBase();
virtual                     ~PdfStreamBase();
        unsigned int        GetPos() { return fNumBytes; }
virtual int                 Write(const void* ptr, int count) = 0;
        PdfStreamBase&      operator<<(const char c)
                                { Write(&c, 1); return *this; }
        PdfStreamBase&      operator<<(unsigned char c)
                                { return (*this) << (char)c; }
        PdfStreamBase&      operator<<(const char* s);
        PdfStreamBase&      operator<<(unsigned int i);
        PdfStreamBase&      operator<<(int i);
        PdfStreamBase&      operator<<(float f)
                                { return (*this) << (double)f; }
        PdfStreamBase&      operator<<(double d);
protected:
        unsigned int        fNumBytes;
        char                fDecimalPoint;
};

/*----------------------------------------------------------------------------*/
/*----- PdfFileStream class --------------------------------------------------*/

class PdfFileStream : public PdfStreamBase
{
public:
                            PdfFileStream(FILE *file, bool closefile = false);
                            PdfFileStream(const char *filename);
virtual                     ~PdfFileStream();
virtual int                 Write(const void* ptr, int count);
private:
        FILE*               fFile;
        bool                fCloseFile;
        int                 fErrorNo;
};

/*----------------------------------------------------------------------------*/
/*----- PdfEncryptor class ---------------------------------------------------*/

#ifdef USE_ENCRYPTION
class PdfEncryptor
{
public:
                        PdfEncryptor(const unsigned char* key);
                        ~PdfEncryptor();
        void            Init(PdfOID id, unsigned int gen_no);
        void            Reset();
        void            CryptBuf(const unsigned char* src,
                            unsigned char* dst, unsigned int len);
private:
        unsigned char   fKey[PDF_ENCRYPT_KEY_LEN + 5];
        unsigned char   fDigest[PDF_MD5_KEY_LEN];
        unsigned char*  fARC4Key;
};
#else
typedef void PdfEncryptor;
#endif

/*----------------------------------------------------------------------------*/
/*----- PdfMemStream class ---------------------------------------------------*/

class PdfMemStream : public PdfStreamBase
{
public:
                        PdfMemStream(const int bufSize = PDF_DEF_BUF_SIZE);
virtual                 ~PdfMemStream();
virtual int             Write(const void* ptr, int count);
        size_t          GetSize() { 
                            return (fBuf == NULL) ? 0 :
                                (fBuf->CountItems() - 1) * fBufSize + 
                                fCurrentPos; 
                        }
        size_t          WriteToStream(PdfStreamBase* out, 
                            PdfEncryptor* eobj = NULL);
        size_t          WriteToStreamDeflate(PdfStreamBase* out,
                            PdfEncryptor* eobj = NULL);
        void            FreeData();
        unsigned int    BufSize()    { return fBufSize; }
        const void*     GetBuf(int index, unsigned int* size);
        int             GetBufCount() {
                            return (fBuf == NULL) ? 0 : (fBuf->CountItems());
                        }
		PdfMemStream*   Duplicate() const;
private:
        void            InternalWrite(unsigned char** ptr, int *count);
        void            InitStream();
        PdfList*        fBuf;
        unsigned int    fBufSize;
        unsigned int    fCurrentPos;
        unsigned char*  fCurrentPtr;
};

/*----------------------------------------------------------------------------*/

typedef enum pdf_object_type_enum
{
    PDF_OBJ_TYPE_DIRECT,
    PDF_OBJ_TYPE_INDIRECT,
    PDF_OBJ_TYPE_VIRTUAL
} pdf_object_type;

typedef enum pdf_object_class_enum
{
    ocUnknown,
    ocBoolean,
    ocNumber,
    ocReal,
    ocText,
    ocBinary,
    ocName,
    ocArray,
    ocDictionary,
    ocStream,
    ocUnicodeText,
    ocNull
} pdf_object_class;

/*----------------------------------------------------------------------------*/
/*----- PdfObject class ------------------------------------------------------*/

class PdfObject
{
public:
                            PdfObject();
                            PdfObject(PdfOID objectID);
virtual                     ~PdfObject();
virtual pdf_object_type     GetObjectType()     { return (fObjectID > 0 ?
                                PDF_OBJ_TYPE_INDIRECT : PDF_OBJ_TYPE_DIRECT); }
virtual pdf_object_class    GetClass()          { return ocUnknown; }
        PdfOID              GetObjectID()       { return fObjectID; }
        unsigned int        GetGenerationNo()   { return fGenerationNo; }
virtual bool                IsLockedObject()    { return fObjectID < 0; }
        void                SetObjectID(PdfOID value);
        void                WriteToStream(PdfStreamBase* out,
                                PdfEncryptor* e)
                                { InternalWriteStream(out, e); }
        void                WriteValueToStream(PdfStreamBase* out, 
                                PdfEncryptor* e);
        size_t              WriteEscapeName(PdfStreamBase* stream, 
                                const char* text);
        size_t              WriteEscapeText(PdfStreamBase* stream, 
                                const char* text);
protected:
        void                WriteBinary(PdfStreamBase* out, unsigned char* buf,
                                unsigned int len, PdfEncryptor* e);
virtual void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e) = 0;
private:
        PdfOID              fObjectID;
        unsigned int        fGenerationNo;
};

/*----------------------------------------------------------------------------*/
/*----- PdfVirtualObject class -----------------------------------------------*/

class PdfVirtualObject : public PdfObject
{
public:
                        PdfVirtualObject(PdfOID objectID) : 
                            PdfObject(objectID) {}
    pdf_object_type     GetObjectType()     { return PDF_OBJ_TYPE_VIRTUAL; }
    bool                IsLockedObject()    { return true; }
protected:
    void                InternalWriteStream(PdfStreamBase* out,
                            PdfEncryptor* e) {
                            *out << GetObjectID() 
                                 << " " 
                                 << GetGenerationNo() 
                                 << " R"; 
                        } 
};

/*----------------------------------------------------------------------------*/
/*----- PdfNullObject class --------------------------------------------------*/

class PdfNullObject : public PdfObject
{
public:
                            PdfNullObject() : PdfObject() {}
                            PdfNullObject(PdfOID objectID) : 
                                PdfObject(objectID) {}
        pdf_object_class    GetClass()          { return ocNull; }
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e)    { *out << "null"; }
};

/*----------------------------------------------------------------------------*/
/*----- PdfBoolean class -----------------------------------------------------*/

class PdfBoolean : public PdfObject
{
public:
                            PdfBoolean() : PdfObject()  { fValue = false; }
                            PdfBoolean(PdfOID objectID, bool value)
                                : PdfObject(objectID)       { fValue = value; }
                            PdfBoolean(bool value) : PdfObject() 
                                { fValue = value; }
        pdf_object_class    GetClass()              { return ocBoolean; }
        bool                GetValue()              { return fValue; }
        void                SetValue(bool value)    { fValue = value; }
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e)
                                { fValue ? *out << "true" : *out << "false"; }
private:
        bool                fValue;
};

/*----------------------------------------------------------------------------*/
/*----- PdfNumber class ------------------------------------------------------*/

class PdfNumber : public PdfObject
{
public:
                            PdfNumber() : PdfObject()       { fValue = 0; }
                            PdfNumber(PdfOID objectID, int value)
                                : PdfObject(objectID)       { fValue = value; }
                            PdfNumber(int value)
                                : PdfObject()       { fValue = value; }
        pdf_object_class    GetClass()              { return ocNumber; }
        int                 GetValue()              { return fValue; }
        void                SetValue(int value)           { fValue = value; }
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e)
                                { *out << fValue; }
private:
        int                 fValue;
};

/*----------------------------------------------------------------------------*/
/*----- PdfReal class --------------------------------------------------------*/

class PdfReal : public PdfObject
{
public:
                            PdfReal() : PdfObject()         { fValue = 0; }
                            PdfReal(PdfOID objectID, double value)
                                : PdfObject(objectID)       { fValue = value; }
                            PdfReal(double value)
                                : PdfObject()       { fValue = value; }
        pdf_object_class    GetClass()              { return(ocReal); }
        double              GetValue()              { return(fValue); }
        void                SetValue(double value)  { fValue = value; }
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e)
                                { *out << fValue; }
private:
        double              fValue;
};

/*----------------------------------------------------------------------------*/
/*----- PdfName class --------------------------------------------------------*/

class PdfName : public PdfObject
{
public:
                            PdfName() : PdfObject()     { fValue[0] = 0x00; }
                            PdfName(PdfOID objectID, const char* value);
                            PdfName(const char* value);
        pdf_object_class    GetClass()          { return(ocName); }
const   char*               GetValue()          { return fValue; }
        void                SetValue(const char* value);
        bool                EqualTo(const char* value)
                                { return (strcmp(fValue, value) == 0); }
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e)
                                { WriteEscapeName(out, fValue); }
private:
        char                fValue[PDF_LIMIT_MAX_NAME + 1];
};

/*----------------------------------------------------------------------------*/
/*----- PdfText class --------------------------------------------------------*/

class PdfText : public PdfObject
{
public:
                            PdfText() : PdfObject()     { fValue = NULL; }
                            PdfText(PdfOID objectID, const char* value);
                            PdfText(const char* value);
                            ~PdfText();
        pdf_object_class    GetClass()           { return(ocText); }
const   char*               GetValue()           { return fValue; }
        void                SetValue(const char* value);
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e);
private:
        char*               fValue;
};

/*----------------------------------------------------------------------------*/
/*----- PdfBinary class ------------------------------------------------------*/

class PdfBinary : public PdfObject
{
public:
                            PdfBinary() : PdfObject() { fData = NULL; }
                            PdfBinary(PdfOID objectID)
                                : PdfObject(objectID) { fData = NULL; };
        pdf_object_class    GetClass()          { return ocBinary; }
        void                SetData(const void* data, unsigned int length,
                                bool encryptable = true);
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e);
private:
        const void*         fData;
        unsigned int        fLength;
        bool                fEncryptable;
};

/*----------------------------------------------------------------------------*/
/*----- PdfXrefEntry class ---------------------------------------------------*/

class PdfXrefEntry
{
public:
                    PdfXrefEntry(PdfObject *object);
                    ~PdfXrefEntry();
    PdfEntryType    GetEntryType()      { return fEntryType; };
    int             GetByteOffset()     { return fByteOffset; };
    unsigned int    GetGenerationNo()   { return fGenerationNo; };
    void            SetByteOffset(int byteOffset) { fByteOffset = byteOffset; };
    void            SetGenerationNo(unsigned int generationNo)
                        { fGenerationNo = generationNo; };
    void            SetEntryType(PdfEntryType entryType)   
                        { fEntryType = entryType; };
    void            SetObject(PdfObject* object);
    bool            HasObject()         { return (fObject != NULL); };
    PdfObject*      GetObject()         { return fObject; };
private:
    PdfEntryType    fEntryType;
    int             fByteOffset;
    unsigned int    fGenerationNo;
    PdfObject*      fObject;
};

/*----------------------------------------------------------------------------*/
/*----- PdfXref class --------------------------------------------------------*/

class PdfDoc;

class PdfXref 
{
public:
                    PdfXref(PdfDoc* doc);
                    ~PdfXref();
    void            Clear();
    int             GetCount()      { return fEntries->CountItems(); };
    PdfObject*      GetObject(PdfOID objectID);
    PdfOID          AddObject(PdfObject* object);
    void            WriteToStream(PdfStreamBase* out, PdfEncryptor* e);
    int             GetAddr()       { return fAddr; };
    PdfDoc*         GetDoc()        { return fDoc; };
    void            SetError(int err);
protected:
    PdfXrefEntry*   GetEntry(int index) {
                        return (fEntries == NULL) ? NULL :
                            (PdfXrefEntry*)fEntries->ItemAt(index); 
                    };
private:
    void            Init();
    PdfList*        fEntries;
    int             fAddr;
    PdfDoc*         fDoc;
};

/*----------------------------------------------------------------------------*/
/*----- PdfArray class -------------------------------------------------------*/

class PdfArray : public PdfObject
{
public:
                            PdfArray(PdfXref* xref);
                            PdfArray(PdfOID objectID, PdfXref* xref);
                            ~PdfArray();
        pdf_object_class    GetClass()          { return(ocArray); }
        int                 GetCount()  { 
                                return (fItems == NULL) ? 0 : 
                                fItems->CountItems(); 
                            }
        PdfObject*          GetItem(int index);
        int                 GetAsInteger(int index);
        double              GetAsReal(int index);
        void                Add(PdfObject* value);
        void                Add(const pdf_box value);
        void                Add(const pdf_rect value);
        void                Insert(PdfObject* value, int index);
        int                 IndexOf(PdfObject* value) { 
                                return (fItems == NULL) ? 0 : 
                                fItems->IndexOf(value); 
                            }
        PdfXref*            GetXref()      { return(fXref); }
        void                Clear();
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e);
private:
        void                CheckList() { 
                                if (fItems == NULL) 
                                    fItems = new PdfList(); 
                            }
        PdfList*            fItems;
        PdfXref*            fXref;
};

/*----------------------------------------------------------------------------*/
/*----- PdfDictElement class -------------------------------------------------*/

class PdfDictElement
{
public:
                    PdfDictElement();
                    PdfDictElement(PdfName* key, PdfObject* value);
                    PdfDictElement(const char* key, PdfObject* value);
                    ~PdfDictElement();
    void            SetValue(PdfObject* value);
    PdfObject*      GetValue()              { return(fValue); }
    PdfName*        GetKey()                { return(fKey); }
private:
    PdfName*        fKey;
    PdfObject*      fValue;
};

/*----------------------------------------------------------------------------*/
/*----- PdfDictionary class --------------------------------------------------*/

class PdfDictionary : public PdfObject
{
public:
                            PdfDictionary(PdfXref* xref);
                            PdfDictionary(PdfOID objectID, PdfXref* xref);
                            ~PdfDictionary();
        pdf_object_class    GetClass()          { return(ocDictionary); }
        int                 GetCount();
        PdfObject*          GetValue(const char* key);
const   char*               GetKeyValue(unsigned int index);
        void                AddElement(const char* key, PdfObject* value);
        void                RemoveElement(const char* key);
        void                RemoveElement(unsigned int index);
        PdfXref*            GetXref()      { return(fXref); }
        bool                IsTypeOf(const char* type);
const   char*               GetTextValue(const char* key);
const   char*               GetNameValue(const char* key);
        void                SetError(int err);
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e);
        void                SetType(const char* type) {
                                AddElement("Type", new PdfName(type)); 
                            }
private:
        void                CheckList() { 
                                if (fElements == NULL) 
                                    fElements = new PdfList();
                            }
        PdfDictElement*     GetElementAt(int idx) { 
                                return (PdfDictElement*)fElements->ItemAt(idx); 
                            }
        PdfDictElement*     GetElement(const char* key);
        PdfList*            fElements;
        PdfXref*            fXref;
};

/*----------------------------------------------------------------------------*/
/*----- PdfStream class ------------------------------------------------------*/

class PdfStream : public PdfDictionary
{
public:
                            PdfStream(PdfXref* xref);
                            PdfStream(PdfOID objectID, PdfXref* xref);
                            ~PdfStream();
        void                AddFilter(pdf_filter filter);
        void                RemoveFilter(pdf_filter filter);
        pdf_object_class    GetClass()          { return ocStream; }
        PdfMemStream*       GetStream();
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e);
private:
        PdfMemStream*       fStream;
        int                 fFilter;
};

/*----------------------------------------------------------------------------*/
/*----- pdf_char_data struct -------------------------------------------------*/

struct pdf_char_data_struct {
    int char_cd;
    char* char_name;
    int width;
};
typedef pdf_char_data_struct pdf_char_data;

/* read only data */
struct pdf_char_data_struct_ro {
    int char_cd;
    const char* char_name;
    int width;
};
typedef pdf_char_data_struct_ro pdf_char_data_ro;

/*----------------------------------------------------------------------------*/
/*----- PdfAutoPtrObject class -----------------------------------------------*/

class PdfAutoPtrMgr;

class PdfAutoPtrObject
{
friend class PdfAutoPtrMgr;
public:
                        PdfAutoPtrObject();
virtual                 ~PdfAutoPtrObject();
virtual pdf_auto_ptr_object_type    
                        GetType()   { return PDF_OBJECT_UNKNOWN; }
        bool            IsValid()   { return fMgr != NULL; }
protected:
virtual void            Init()  {};
        
private:
        PdfAutoPtrMgr*  fMgr;
};       

class PdfAutoPtrMgr
{
public:
                        PdfAutoPtrMgr();
                        ~PdfAutoPtrMgr();
        void            RegisterObject(PdfAutoPtrObject* obj);
        void            UnRegisterObject(PdfAutoPtrObject* obj);
        unsigned int    CountObjects()  { return fList->CountItems(); }
private:
        PdfList*        fList;
};

/*----------------------------------------------------------------------------*/
/*----- PdfType1FontDef class ------------------------------------------------*/

class PdfType1FontDef : public PdfAutoPtrObject
{
public:
                        	PdfType1FontDef();
                        	PdfType1FontDef(const char* basefont);
virtual                 	~PdfType1FontDef();
        unsigned int    	Widths(const char* char_name);
        unsigned int    	Widths(unsigned char c);

        /* font-specific encording only. */     
        int             	FirstChar()             { return fFirstChar; }
        int             	LastChar()              { return fLastChar; } 

        const char*     	BaseFont()              { return fBaseFont; }
        int             	Ascent()                { return fAscent; }        
        int             	Descent()               { return fDescent; }
        unsigned int    	Flags()                 { return fFlags; }
        const char*     	FontName()              { return fFontName; }
        pdf_box         	FontBBox()              { return fFontBBox; }
        int             	ItalicAngle()           { return fItalicAngle; }
        unsigned int    	StemV()                 { return fStemV; }
        int             	AvgWidth()              { return fAvgWidth; }
        int             	Leading()               { return fLeading; }
        int             	MaxWidth()              { return fMaxWidth; }
        int             	MissingWidth()          { return fMissingWidth; }
        int             	StemH()                 { return fStemH; }
        int             	XHeight()               { return fXHeight; }
        int             	CapHeight()             { return fCapHeight; }
        const char*     	CharSet()               { return fCharSet; }
        bool            	IsBase14Font()          { return fIsBase14Font; }
        unsigned int    	Length1()               { return fLength1; }
        unsigned int    	Length2()               { return fLength2; }
        unsigned int    	Length3()               { return fLength3; }
        PdfMemStream*   	FontData()              { return fFontData; }
        PdfDictionary*  	Descriptor()            { return fDescriptor; }
        void            	SetDescriptor(PdfDictionary* desc)
                            	        { fDescriptor = desc; }
        pdf_encoding    	DefaultEncoding()   
                            	        { return fDefaultEncoding; }
        void            	Clear();
        void            	LoadFromFile(const char* afmfile, 
								const char* fontfile);
		PdfType1FontDef* 	Duplicate() const;
protected:
        void            SetWidths(const pdf_char_data_ro* array);
        void            FreeWidths();
        char            fBaseFont[PDF_LIMIT_MAX_NAME + 1];
        pdf_char_data*  fWidths;
        int             fWidthsCount;
        int             fFirstChar;
        int             fLastChar;
        int             fAscent;
        int             fDescent;
        unsigned int    fFlags;
        char            fFontName[PDF_LIMIT_MAX_NAME + 1];
        pdf_box         fFontBBox; 
        int             fItalicAngle;
        unsigned int    fStemV;
        int             fAvgWidth;
        int             fLeading;
        int             fMaxWidth;
        int             fMissingWidth;
        int             fStemH;
        int             fXHeight;
        int             fCapHeight;
        char            fCharSet[PDF_LIMIT_MAX_NAME + 1];
        unsigned int    fLength1;
        unsigned int    fLength2;
        unsigned int    fLength3;
        bool            fIsBase14Font;
        pdf_encoding    fDefaultEncoding;
        PdfMemStream*   fFontData;
        PdfDictionary*  fDescriptor;
private:
        void            InitAttributes();
        void            SetParam(char* param, const char* text);
        void            LoadFontFile(const char* filename);
        void            LoadAfmFile(const char* filename);
        char*           GetStrParam(char* str, char* param, int len);
        char*           GetIntParam(char* str, int* param);
        int             FindBuf(const char* findbuf, const char* srcbuf, 
                            int len);
};

/*----------------------------------------------------------------------------*/
/*----- Helvetica Font -------------------------------------------------------*/

class PdfHelveticaFontDef : public PdfType1FontDef
{
public:
                    PdfHelveticaFontDef();
};

class PdfHelveticaBoldFontDef : public PdfType1FontDef
{
public:
                    PdfHelveticaBoldFontDef();
};

class PdfHelveticaObliqueFontDef : public PdfType1FontDef
{
public:
                    PdfHelveticaObliqueFontDef();
};

class PdfHelveticaBoldObliqueFontDef : public PdfType1FontDef
{
public:
                    PdfHelveticaBoldObliqueFontDef();
};

/*----------------------------------------------------------------------------*/
/*----- Times Font -----------------------------------------------------------*/

class PdfTimesRomanFontDef : public PdfType1FontDef
{
public:
                    PdfTimesRomanFontDef();
};

class PdfTimesBoldFontDef : public PdfType1FontDef
{
public:
                    PdfTimesBoldFontDef();
};

class PdfTimesItalicFontDef : public PdfType1FontDef
{
public:
                    PdfTimesItalicFontDef();
};

class PdfTimesBoldItalicFontDef : public PdfType1FontDef
{
public:
                    PdfTimesBoldItalicFontDef();
};

/*----------------------------------------------------------------------------*/
/*----- Courier Font ---------------------------------------------------------*/

class PdfCourierFontDef : public PdfType1FontDef
{
public:
                    PdfCourierFontDef();
};

class PdfCourierBoldFontDef : public PdfType1FontDef
{
public:
                    PdfCourierBoldFontDef();
};

class PdfCourierObliqueFontDef : public PdfType1FontDef
{
public:
                    PdfCourierObliqueFontDef();
};

class PdfCourierBoldObliqueFontDef : public PdfType1FontDef
{
public:
                    PdfCourierBoldObliqueFontDef();
};

/*----------------------------------------------------------------------------*/
/*----- Symbol Font ----------------------------------------------------------*/

class PdfSymbolFontDef : public PdfType1FontDef
{
public:
                    PdfSymbolFontDef();
};

/*----------------------------------------------------------------------------*/
/*----- ZapfDingbats Font ----------------------------------------------------*/

class PdfZapfDingbatsFontDef : public PdfType1FontDef
{
public:
                    PdfZapfDingbatsFontDef();
};

/*----------------------------------------------------------------------------*/
/*------ PdfEncodingDef class ------------------------------------------------*/

struct pdf_differences_element {
    unsigned int code;
    const char* char_name;
};

class PdfEncodingDef : public PdfAutoPtrObject
{
public:
                        PdfEncodingDef();
virtual                 ~PdfEncodingDef();
        int             GetCharCode(unsigned short unicode);
        int             GetCharCode(const char* char_name);
        const char*     GetCharName(unsigned int code);
        unsigned short  GetUnicode(unsigned int code);
        pdf_encoding    BaseEncoding()      { return fBaseEncoding; }
virtual PdfObject*      GetEncoding(PdfXref* xref);
        unsigned int    FirstChar()         { return fFirstChar; }
        unsigned int    LastChar()          { return fLastChar; }
        void            OverrideUnicodeArray(const unsigned short* array);
        int             ToUnicode(const char* src, unsigned char* dst,
                            int* len);
static  pdf_encoding    StringToPdfEncoding(const char* encoding);
protected:
        void            AddDifference(unsigned int code, const char* char_name);
        void            SetUnicodeArray(unsigned int first_char,
                            unsigned int last_char, 
                            const unsigned short* array);
        const char*     PdfEncodingToString(pdf_encoding encoding);
        const char*     UnicodeToGryphName(unsigned short unicode);
        unsigned short  GryphNameToUnicode(const char* gryph_name);
        unsigned int    fFirstChar;
        unsigned int    fLastChar;
        pdf_encoding    fBaseEncoding;
private:
        void            ClearDefferences();
        PdfList*        fDifferences;
        unsigned short  fUnicodeArray[255];
};

/*----------------------------------------------------------------------------*/
/*------ predefined encoding -------------------------------------------------*/

class PdfStandardEncoding : public PdfEncodingDef
{
public:
                        PdfStandardEncoding();
};

class PdfWinAnsiEncoding : public PdfEncodingDef
{
public:
                        PdfWinAnsiEncoding();
};

class PdfMacRomanEncoding : public PdfEncodingDef
{
public:
                        PdfMacRomanEncoding(); 
};

class PdfSymbolFontEncoding : public PdfEncodingDef
{
public:
                        PdfSymbolFontEncoding();
        PdfObject*      GetEncoding(PdfXref* xref)      { return NULL; }
};

class PdfZapfDingbatsFontEncoding : public PdfEncodingDef
{
public:
                        PdfZapfDingbatsFontEncoding();
        PdfObject*      GetEncoding(PdfXref* xref)      { return NULL; }
};


/*----------------------------------------------------------------------------*/
/*----- PdfFontBase class ----------------------------------------------------*/

class PdfFontBase : public PdfStream
{
public:
                            PdfFontBase(PdfXref* xref);
                            ~PdfFontBase();
        const char*         Name()    { return fName; }
virtual pdf_text_width      TextWidth(const char* text);
virtual pdf_text_width      TextWidths(const char* text, 
                                unsigned int* widths);
virtual unsigned int        MeasureText(const char* text, double width,
                                double fontsize, double charspace,
                                double wordspace, double* realwdth = NULL);
virtual pdf_writing_mode    WritingMode()   { return fWritingMode; }
virtual int                 Ascent() = 0;
virtual int                 Descent() = 0;
protected:
        void                SetWritingMode(pdf_writing_mode mode)
                                { fWritingMode = mode; }
        char*               fName;
        pdf_writing_mode    fWritingMode;
private:
};

typedef PdfFontBase PdfFont;

/*----------------------------------------------------------------------------*/
/*----- PdfFontMgr class -----------------------------------------------------*/

class PdfFontMgr
{
public:
                        PdfFontMgr(PdfXref *xref);
                        ~PdfFontMgr();
        int             RegisterFont(PdfFont* font);
        PdfFont*        GetFont(const char* name);
        PdfFont*        GetFont(unsigned int index)
                            { return (PdfFont*)fList->ItemAt(index); }
        unsigned int    CountFonts()  { return fList->CountItems(); }
private:
        PdfList*        fList;
        PdfXref*        fXref;
};

/*----------------------------------------------------------------------------*/
/*----- PdfType1Font class ---------------------------------------------------*/

class PdfType1Font : public PdfFont
{
public:
                            PdfType1Font(PdfXref* xref);
                            ~PdfType1Font();
        pdf_text_width      TextWidth(const char* text);
        pdf_text_width      TextWidths(const char* text, unsigned int* widths);
        unsigned int        MeasureText(const char* text, double width,
                                double fontsize, double charspace, 
                                double wordspace, double* realwidth = NULL);
        unsigned int        CharWidth(const unsigned char c)  {
                                return (fFirstChar <= c && c <= fLastChar) ? 
                                    fWidths[c - fFirstChar] : fMissingWidth; }
        PdfEncodingDef*     Encording()     { return fEncoding; }
        PdfType1FontDef*    FontDef()       { return fFontDef; }
        bool                IsBase14Font()  { return fFontDef->IsBase14Font(); }
        bool                GetValid()      { return fValid; }
        void                SetAttributes(const char* name,
                                PdfType1FontDef* fontdef,
                                PdfEncodingDef* encoding);
        int                 Ascent();
        int                 Descent();
protected:      
        bool                CreateDescriptor();
        bool                fValid;
private:
        PdfDictionary*      fDescriptor;
        PdfEncodingDef*     fEncoding;
        PdfType1FontDef*    fFontDef;
        int                 fFirstChar;
        int                 fLastChar;
        int                 fMissingWidth;
        int*                fWidths;
		int					fFlags;
};

/*---------------------------------------------------------------------------*/
/*----- PdfCMap -------------------------------------------------------------*/

class PdfCIDType2FontDef;
typedef PdfCIDType2FontDef PdfCIDFontDef;

class PdfCIDType2Font;
typedef PdfCIDType2Font PdfCIDFont;

class PdfCMap : public PdfAutoPtrObject
{
public:
                            PdfCMap();
virtual                     ~PdfCMap();
        pdf_cid             GetCID(unsigned int code);
        unsigned int        GetUnicode(unsigned int code);
virtual const char*         GetCMapName() = 0;
virtual void                ParseText(const char* text, 
                                pdf_byte_type* btype) = 0;
        int                 ToUnicode(const char* src, unsigned char* dst, 
                                int* len);
        void                AddCIDSystemInfo(PdfCIDFont* font);
virtual pdf_writing_mode    GetWritingMode()    
                                { return PDF_WMODE_HORIZONTAL; }
protected:
        void                AddCMap(const pdf_cid_range* range);
        void                SetCIDSystemInfo(const char* registry, 
                                const char* ordering, unsigned int supplement);
        void                SetUnicodeArray(const pdf_mb_unicode_map1* array1,
                                const pdf_mb_unicode_map2* array2);
private:
    unsigned short          fUnicodeArray[256][256];
    pdf_cid                 fCMapArray[256][256];
    char*                   fRegistry;
    char*                   fOrdering;
    unsigned int            fSupplement;
};

/*----------------------------------------------------------------------------*/
/*----- PdfUnicodeText -------------------------------------------------------*/

class PdfUnicodeText : public PdfBinary
{
public:
                            PdfUnicodeText();
                            PdfUnicodeText(PdfCMap* cmap);
                            PdfUnicodeText(PdfEncodingDef* def);
                            ~PdfUnicodeText();
        pdf_object_class    GetClass()           { return(ocUnicodeText); }
        void                SetText(const char* text);
        const char*         GetText();
private:
        unsigned char*      fValue;
        PdfCMap*            fCMap;
        PdfEncodingDef*     fDef;
};

/*---------------------------------------------------------------------------*/
/*----- PdfType0Font --------------------------------------------------------*/

class PdfType0Font : public PdfFont
{
public:
                            PdfType0Font(PdfXref* xref);
                            ~PdfType0Font();
        pdf_text_width      TextWidth(const char* text) {
                                return TextWidths(text, NULL);
                            }
        pdf_text_width      TextWidths(const char* text, unsigned int* widths);
        unsigned int        MeasureText(const char* text, double width,
                                double fontsize, double charspace,
                                double wordspace, double* realwidth = NULL);
        PdfCIDType2Font*    DescendantFont()    { return fDescendantFont; }
        void                SetAttributes(const char* name, PdfCIDFont* font, 
                                PdfCMap* cmap);
        bool                GetValid()      { return fValid; }
        int                 Ascent();
        int                 Descent();
private:
        PdfCIDType2Font*    fDescendantFont;
        PdfCMap*            fCMap;
        bool                fValid;
};

/*---------------------------------------------------------------------------*/
/*----- PdfCIDType2FontDef --------------------------------------------------*/

class PdfCIDType2FontDef : public PdfAutoPtrObject
{
friend class PdfCIDType2Font;
public:
                            PdfCIDType2FontDef();
virtual                     ~PdfCIDType2FontDef();
        unsigned int        CIDWidth(pdf_cid cid);
        const char*         BaseFont()              { return fBaseFont; }
        int                 Ascent()                { return fAscent; }
        int                 Descent()               { return fDescent; }
        int                 CapHeight()             { return fCapHeight; }
        unsigned int        Flags()                 { return fFlags; }
        pdf_box             FontBBox()              { return fFontBBox; }
        int                 ItalicAngle()           { return fItalicAngle; }
        unsigned int        StemV()                 { return fStemV; }
        int                 DW()                    { return fDW; }
        int*                DW2()                   { return fDW2; }
        int                 MissingWidth()          { return fMissingWidth; }
        pdf_cid_width*      GetWidths(int index) {
                                return (fWidths == NULL) ? 0 :
                                    (pdf_cid_width*)fWidths->ItemAt(index);
                            }
        int                 NumWidths();
protected:
        void                AddWidths1(pdf_cid fromcid, pdf_cid tocid,
                                unsigned int width);
        void                AddWidths2(pdf_cid fromcid, pdf_cid tocid,
                                const unsigned int* widths);
        void                SetBaseFont(const char* value) 
                                { SetParam(&fBaseFont, value); }
virtual void                Init() = 0;
        int                 fAscent;
        int                 fDescent;
        int                 fCapHeight;
        unsigned int        fFlags;
        pdf_box             fFontBBox;
        int                 fItalicAngle;
        unsigned int        fStemV;
        int                 fAvgWidth;
        int                 fLeading;
        int                 fMaxWidth;
        int                 fMissingWidth;
        int                 fStemH;
        int                 fDW;
        int                 fDW2[2];
private:
        void                SetParam(char** dst, const char* src);
        PdfList*            fWidths;
        char*               fBaseFont;
};

/*---------------------------------------------------------------------------*/
/*----- PdfCIDType2Font -----------------------------------------------------*/

class PdfCIDType2Font : public PdfFont
{
public:
                            PdfCIDType2Font(PdfXref* xref);
                            ~PdfCIDType2Font();
        PdfCIDFontDef*      FontDef()               { return fFontDef; }
        const char*         BaseFont() {
                                return fValid ? fFontDef->BaseFont() : NULL;
                            }
        bool                GetValid()              { return fValid; }
        unsigned int        CIDWidth(pdf_cid cid) {
                                return fValid ? fFontDef->CIDWidth(cid) : 0;
                            }
        void                SetAttributes(PdfCIDFontDef* fontdef);
        int                 Ascent();
        int                 Descent();
private:
        void                CreateDescriptor();
        bool                fValid;
        PdfCIDFontDef*      fFontDef;
        PdfDictionary*      fDescriptor;
};
        
/*----------------------------------------------------------------------------*/
/*----- PdfPageBase class ----------------------------------------------------*/

enum pdf_pages_type {
    ptHasNoKids = 0,
    ptHasPages,
    ptHasPage
};

class PdfPages;
class PdfPage;
class PdfXObjectBase;
typedef PdfXObjectBase PdfXObject;

class PdfPageBase : public PdfDictionary
{
    friend class PdfPages;
public:
                        PdfPageBase(PdfXref* xref);
                        ~PdfPageBase();
virtual void            Init();
        pdf_box         MediaBox();
        pdf_box         CropBox();
        PdfDictionary*  Resources();
        int             Rotate();
        PdfPages*       Parent()        { return fParent; }
        PdfDictionary*  GetResource(const char* element_name);
        void            SetMediaBox(pdf_box rect);
        void            SetCropBox(pdf_box rect);
        void            SetResources(PdfDictionary* res);
        void            SetRotate(int rotate);
        int             CountFonts();
        bool            AddFont(PdfFont* font, const char* fontname);
        const char*     GetFontName(PdfFont* font);
        bool            AddXObject(PdfXObject* xobject, const char* name);
        const char*     GetXObjectName(PdfXObject* xobject);
        int             GetIndex();
        void            SetSize(int width, int height);
        void            AddProcSet(int procset) {
                            if ((fProcSet & procset) != procset)
                                InternalAddProcSet(procset); 
                        } 
protected:
        PdfObject*      FindElement(const char* name);
        pdf_box         GetElementRect(const char* name);
        void            SetElementRect(const char* name, pdf_box rect);
virtual unsigned int    GetPageCount() = 0;
        PdfPages*       fParent;
private:
        void            InternalAddProcSet(int procset);
        int             fProcSet;
};

/*----------------------------------------------------------------------------*/
/*----- PdfPages class -------------------------------------------------------*/

class PdfXObjectMgr;

class PdfPages : public PdfPageBase
{
public:
                        PdfPages(PdfXref *xref, PdfFontMgr* fmgr = NULL,
                            PdfXObjectMgr* omgr = NULL);
        void            Init();
        PdfFontMgr*     FontMgr()       { return fFontMgr; }
        PdfXObjectMgr*  XObjectMgr()    { return fXObjectMgr; }
        PdfPage*        AddPage(int index = -1);
        PdfPages*       AddPages(int index = -1);
        unsigned int    GetKidsCount();
        int             IndexOf(PdfPageBase* page);
protected:
        void            InternalWriteStream(PdfStreamBase* out,
                            PdfEncryptor* e);
private:
        void            SetParent(PdfPages *parent);
        PdfArray*       Kids()      { return (PdfArray*)GetValue("Kids"); }
        bool            AddKids(PdfPageBase* page);
        bool            InsertKids(PdfPageBase* page, int index);
        unsigned int    GetPageCount();
        PdfFontMgr*     fFontMgr;
        PdfXObjectMgr*  fXObjectMgr;
        pdf_pages_type  fPagesType;
};

/*----------------------------------------------------------------------------*/
/*----- PdfPage class --------------------------------------------------------*/

class PdfContents;
class PdfAnnotation;
class PdfLinkAnnot;
class PdfTextAnnot;
class PdfDestination;

class PdfPage : public PdfPageBase
{
    friend class PdfPages;
public:
                    PdfPage(PdfXref* xref)
                        : PdfPageBase(xref)  { fCanvas = NULL; }
    PdfContents*    Canvas();
    PdfFontMgr*     FontMgr()       { return ((PdfPages*)fParent)->FontMgr(); }
    PdfXObjectMgr*  XObjectMgr()    {
                         return ((PdfPages*)fParent)->XObjectMgr(); 
                    }
    int             Width()         { return MediaBox().right; }
    int             Height()        { return MediaBox().top; }
    PdfLinkAnnot*   AddLink(pdf_rect rect, PdfDestination* dest, 
                        pdf_annot_hl_mode mode = PDF_ANNOT_HL_EOF);
    PdfTextAnnot*   AddTextAnnot(pdf_rect rect);
protected:
    void            AddAnnotation(PdfAnnotation* annot);
private:
    void            SetParent(PdfPages *parent);
    unsigned int    GetPageCount()  { return 1; }
    PdfContents*    fCanvas;
};

/*----------------------------------------------------------------------------*/
/*----- PdfDestination class -------------------------------------------------*/

class PdfDestination : public PdfArray
{
public:
                            PdfDestination(PdfPage* page, 
                                    bool create_as_shared = true); 
    PdfPage*                Page()  { return fPage; }
    pdf_destination_type    Type();
    void                    SetXYZ(double left, double top, double zoom);
    void                    SetFit();
    void                    SetFitH(double top);
    void                    SetFitV(double left);
    void                    SetFitR(double left, double bottom, double right, 
                                double top);
    void                    SetFitB();
    void                    SetFitBH(double top);
    void                    SetFitBV(double left);
    double                   GetParam(int index);
    bool                    CheckValid();
    bool                    IsShared()  { return fShared; }
private:
    void                    CheckShared();
    PdfPage*                fPage;
    bool                    fShared;
};

/*---------------------------------------------------------------------------*/
/*----- PdfBorderStyle ------------------------------------------------------*/

class PdfBorderStyle : public PdfDictionary
{
public:
                            PdfBorderStyle(PdfXref *xref);
        pdf_bs_subtype      Subtype()   { return fSubtype; }
        void                SetSubtype(pdf_bs_subtype type) 
                                { fSubtype = type; }
        double               Width()     { return fWidth; }
        void                SetWidth(double width)   
                                { (width >= 0) ? fWidth = width: fWidth = 0; }
        void                SetDash(unsigned int on, unsigned int off,
                                unsigned int phase);
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e);
private:
        pdf_bs_subtype      fSubtype;
        double               fWidth;
        unsigned int        fDashOn;
        unsigned int        fDashOff;
        unsigned int        fPhase; 
};
        
/*---------------------------------------------------------------------------*/
/*----- PdfAnnotation -------------------------------------------------------*/

class PdfAnnotation : public PdfDictionary
{
public:
                            PdfAnnotation(PdfXref *xref, pdf_annot_type type);
        pdf_annot_type      Subtype()   { return fSubtype; }
        pdf_rect            Rect()      { return fRect; }
        void                SetRect(pdf_rect rect); 
        const char*         SubtypeText(pdf_annot_type type);
        void                SetBorder(double horiz, double vert, double width);
        void                GetBorder(double* horiz, double* vert, 
                                double* width);
        PdfBorderStyle*     GetBorderStyle(); 
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e);
private:
        pdf_annot_type      fSubtype;
        pdf_rect            fRect;
        double              fBorder[3];
};

/*----------------------------------------------------------------------------*/
/*----- PdfLinkAnnot ---------------------------------------------------------*/

class PdfLinkAnnot : public PdfAnnotation
{
public:
                        PdfLinkAnnot(PdfXref* xref) 
                            : PdfAnnotation(xref, PDF_ANNOT_LINK) {}
    pdf_annot_hl_mode   HightlightMode();
    PdfDestination*     Dest();
    void                SetDest(PdfDestination* dest);
    void                SetHightlightMode(pdf_annot_hl_mode mode);
};  

/*----------------------------------------------------------------------------*/
/*----- PdfTextAnnot ---------------------------------------------------------*/

class PdfTextAnnot : public PdfAnnotation
{
public:
                            PdfTextAnnot(PdfXref* xref);
    const char*             GetContents();
    void                    SetContents(const char* text, 
                                PdfEncodingDef* encoding = NULL);
    void                    SetContentsMb(const char* text, PdfCMap* cmap);
    bool                    GetOpened()             { return fOpened; }
    void                    SetOpened(bool value)   { fOpened = value; }
    void                    SetIcon(pdf_annot_icon_names icon) 
                                { fIcon = icon; }
    pdf_annot_icon_names    GetIcon()   { return fIcon; }
protected:
    void                    InternalWriteStream(PdfStreamBase* out,
                                    PdfEncryptor* e);
private:
    PdfObject*              fContents;
    bool                    fOpened;
    pdf_annot_icon_names    fIcon;
};

/*----------------------------------------------------------------------------*/
/*----- operator type --------------------------------------------------------*/

typedef enum pdf_ope_type_enum {
    otUnknown = 0,
    otGeneralGState,
    otSpecialGState,
    otPathConstruction,
    otColor,
    otTextState,
    otTextShowing,
    otTextPositioning,
    otMarkdContent,
    otID
} pdf_ope_type;

/*----------------------------------------------------------------------------*/
/*----- PdfContents class ----------------------------------------------------*/

class PdfContents : public PdfStream
{
public:
                        PdfContents(PdfPage *page);
    void                Init();
    PdfPage*            Page()              { return fPage; }
    int                 Width()             { return fPage->Width(); }
    int                 Height()            { return fPage->Height(); }
    pdf_point           CurPoint()          { return fCurPoint; }
    pdf_point           TextPos()           { return fTextPoint; }
    double              TextLeading()       { return fTextLeading; }
    double              FontSize()          { return fFontSize; }
    const char*         FontName()
                            { return ((fFont != NULL) ? fFont->Name() : NULL); }
    PdfFont*            Font()              { return fFont; }
    double              CharSpace()         { return fCharSpace; }
    double              WordSpace()         { return fWordSpace; }
    double              TextWidth(const char* text, int* numchars = NULL,
                            int* numwords = NULL);
    void                CharWidths(const char* chars, double* widths);
    unsigned int        MeasureText(const char* text, double width, 
                            double* realwidth = NULL);
    void                TextRect(const char* text, pdf_rect rect, 
                            int max_len = PDF_DEF_BUF_SIZE);
    void                TextOut(double x, double y, const char* text);
    pdf_rgb_color       RGBFill()           { return fRGBFill; }
    pdf_rgb_color       RGBStroke()         { return fRGBStroke; }
    double              GrayFill()          { return fGrayFill; }
    double              GrayStroke()        { return fGrayStroke; }
    pdf_graphics_mode   GMode()             { return fGMode; }

    /*--- General graphics state ---------------------------------------------*/

    void        SetLineWidth(double linewidth);                          /* w */
    void        SetLineCap(pdf_line_cap_style linecap);                  /* J */
    void        SetLineJoin(pdf_line_join_style linejoin);               /* j */
    void        SetMiterLimit(double miterlimit);                        /* M */
    void        SetDash(unsigned int on, unsigned int off, 
                    unsigned int phase);                                 /* d */
                                                  /* ri --not implemented yet */
    void        SetFlat(unsigned int flatness);                          /* i */
                                                  /* gs --not implemented yet */

    /*--- Special graphic state ----------------------------------------------*/

    void        GSave();                                                 /* q */
    void        GRestore();                                              /* Q */
    void        Concat(double a, double b, double c,
                       double d, double e, double f);                   /* cm */

    /*--- Path construction --------------------------------------------------*/

    void        MoveTo(double x, double y);                              /* m */
    void        LineTo(double x, double y);                              /* l */
    void        CurveTo(double x1, double y1, double x2, 
                    double y2, double x3, double y3);                    /* c */
    void        CurveTo2(double x2, double y2, double x3, double y3);    /* v */
    void        CurveTo3(double x1, double y1, double x3, double y3);    /* y */
    void        ClosePath();                                             /* h */
    void        Rectangle(double x, double y, 
                    double width, double height);                       /* re */

    /*--- Path painting ------------------------------------------------------*/

    void        Stroke();                                                /* S */
    void        ClosePathStroke();                                       /* s */
    void        Fill();                                                  /* f */
    void        Eofill();                                               /* f* */
    void        FillStroke();                                            /* B */
    void        EofillStroke();                                         /* B* */
    void        ClosePathFillStroke();                                   /* b */
    void        ClosePathEofillStroke();                                /* b* */
    void        EndPath();                                               /* n */

    /*--- Clipping paths -----------------------------------------------------*/

    void        Clip();                                                  /* W */
    void        EoClip();                                               /* W* */

    /*--- Text object --------------------------------------------------------*/

    void        BeginText();                                            /* BT */
    void        EndText();                                              /* ET */

    /*--- Text state ---------------------------------------------------------*/

    void        SetCharSpace(double value);                             /* Tc */
    void        SetWordSpace(double value);                             /* Tw */
    void        SetHorizontalScalling(double value);                    /* Tz */
    void        SetTextLeading(double value);                           /* TL */
    void        SetFontAndSize(const char* fontname, double size);      /* Tf */
    void        SetFontAndSize(PdfFont* font, double size);          
    void        SetTextRenderingMode(pdf_text_rendering_mode mode);     /* Tr */
    void        SetTextRaise(double value);                             /* Ts */

    /*--- Text positioning ---------------------------------------------------*/

    void        MoveTextPos(double tx, double ty);                      /* Td */
    void        MoveTextPos2(double tx, double ty);                     /* TD */
    void        SetTextMatrix(double a, double b, double c, 
                    double d, double x, double y);                      /* Tm */
    void        MoveToNextLine();                                       /* T* */

    /*--- Text showing -------------------------------------------------------*/

    void        ShowText(const char* text);                             /* Tj */
                                                                        /* TJ */
    void        ShowTextNextLine(const char* text);                      /* ' */
    void        ShowTextNextLine(double aw, double ac,
                    const char* text);                                   /* " */

    /*--- Color showing ------------------------------------------------------*/

                                                  /* cs --not implemented yet */
                                                  /* CS --not implemented yet */
                                                  /* sc --not implemented yet */
                                                 /* scn --not implemented yet */
                                                  /* SC --not implemented yet */
                                                 /* SCN --not implemented yet */
    void        SetGrayFill(double gray);                                /* g */
    void        SetGrayStroke(double gray);                              /* G */
    void        SetRGBFill(double r, double g, double b);               /* rg */
    void        SetRGBFill(int r, int g, int b);
    void        SetRGBFill(pdf_rgb_color c) 
                    { SetRGBFill(c.red, c.green, c.blue); }
    void        SetRGBStroke(double r, double g, double b);             /* RG */
    void        SetRGBStroke(pdf_rgb_color c)
                    { SetRGBStroke(c.red, c.green, c.blue); }
    void        SetRGBStroke(int r, int g, int b);
    void        SetCMYKFill(double c, double m, double y, double k);     /* k */
    void        SetCMYKStroke(double c, double m, double y, double k);   /* K */

    /*--- Shading patterns ---------------------------------------------------*/
    
                                                  /* sh --not implemented yet */
    
    /*--- In-line images -----------------------------------------------------*/
    
                                                  /* BI --not implemented yet */
                                                  /* ID --not implemented yet */
                                                  /* EI --not implemented yet */

    /*--- XObjects -----------------------------------------------------------*/

    void        ExecuteXObject(const char* name);                       /* Do */
    void        ExecuteXObject(PdfXObject* xobject);

    /*--- Marked content -----------------------------------------------------*/

                                                 /* BMC --not implemented yet */
                                                 /* BDC --not implemented yet */
                                                 /* EMC --not implemented yet */
                                                  /* MP --not implemented yet */
                                                  /* DP --not implemented yet */

    /*--- Compatibility ------------------------------------------------------*/

                                                  /* BX --not implemented yet */
                                                  /* EX --not implemented yet */
protected:
private:
    PdfPage*                fPage;
    PdfFontMgr*             fFontMgr;
    PdfXObjectMgr*          fXObjectMgr;

    double                  fCharSpace;
    double                  fWordSpace;
    double                  fFontSize;
    PdfFont*                fFont;
    double                  fHScalling;
    double                  fTextLeading;
    pdf_text_rendering_mode fRenderingMode;
    double                  fTextRaise;

    double                  fLineWidth;
    pdf_line_cap_style      fLineCap;
    pdf_line_join_style     fLineJoin;
    double                  fMiterLimit;
    unsigned int            fDashOn;
    unsigned int            fDashOff;
    unsigned int            fDashPhase;
    unsigned int            fFlatness;

    pdf_point               fStartPoint;
    pdf_point               fCurPoint;
    pdf_point               fTextPoint;
    pdf_text_matrix         fMatrix;
    pdf_graphics_mode       fGMode;

    pdf_rgb_color           fRGBFill;
    pdf_rgb_color           fRGBStroke;
    double                  fGrayFill;
    double                  fGrayStroke;
};

/*----------------------------------------------------------------------------*/
/*----- PdfCatalog class -----------------------------------------------------*/

class PdfOutlineRoot;
class PdfPageLabel;

class PdfCatalog : public PdfDictionary
{
public:
                    PdfCatalog(PdfXref *xref);
                    PdfCatalog(int objectID, PdfXref *xref);
    void            Init();
    PdfDestination* OpenAction()            { return fOpenAction; }
    PdfOutlineRoot* Outlines();
    pdf_page_layout PageLayout()            { return fPageLayout; }
    pdf_page_mode   PageMode()              { return fPageMode; }
    pdf_page_mode   NonFullScreenPageMode() { return fNonFullScreenPageMode; }
    int             ViewerPreferences()     { return fViewerPreferences; }
    bool            IsSetViewerPreference(int preferences)  { return 
                        ((fViewerPreferences & preferences) == preferences); }
    PdfPages*       Pages();
    void            SetPageLayout(pdf_page_layout layout)
                            { fPageLayout = layout; }
    void            SetPageMode(pdf_page_mode mode)
                            { fPageMode = mode; }
    void            SetNonFullScreenPageMode(pdf_page_mode mode)
                            { fNonFullScreenPageMode = mode; }
    void            SetViewerPreferences(int viewerPreferences)
                            { fViewerPreferences = viewerPreferences; }
    void            SetOpenAction(PdfDestination* action);
    void            AddPageLabel(unsigned int pagenum, pdf_page_num_style style,
                        unsigned int firstpage = 1, const char* prefix = NULL); 
    void            ClearPageLabel();
protected:
    void            InternalAddPageLabel(unsigned int pagenum, 
                        PdfPageLabel* label);
    void            InternalWriteStream(PdfStreamBase* out,
                        PdfEncryptor* e);
private:
    pdf_page_layout fPageLayout;
    pdf_page_mode   fPageMode;
    pdf_page_mode   fNonFullScreenPageMode;
    PdfDestination* fOpenAction;
    int             fViewerPreferences;
};

/*---------------------------------------------------------------------------*/
/*----- PdfPageLabel --------------------------------------------------------*/

class  PdfPageLabel : public PdfDictionary
{
public:
                            PdfPageLabel(PdfXref *xref);
                            ~PdfPageLabel();
        pdf_page_num_style  NumberingStyle()    { return fNumberingStyle; }
        const char*         LabelPrefix()       { return fLabelPrefix; }
        unsigned int        FirstPage()         { return fFirstPage; }
        void                SetNumberingStyle(pdf_page_num_style value)
                                                { fNumberingStyle = value; }
        void                SetLabelPrefix(const char* value);
        void                SetFirstPage(unsigned int page) 
                                                { fFirstPage = page; }
protected:
        void                InternalWriteStream(PdfStreamBase* out,
                                PdfEncryptor* e);
private:
        pdf_page_num_style  fNumberingStyle;
        char*               fLabelPrefix;
        unsigned int        fFirstPage;
};

/*----------------------------------------------------------------------------*/
/*------ PdfInfo -------------------------------------------------------------*/

class PdfInfo : public PdfDictionary
{
public:
                    PdfInfo(PdfXref* xref);
                    PdfInfo(int objectID, PdfXref* xref);
    void            Init();
    pdf_date        CreationDate();
    pdf_date        ModDate();
    const char*     Author()            { return GetTextValue("Author"); }
    const char*     Creator()           { return GetTextValue("Creator"); }
    const char*     Producer()          { return GetTextValue("Producer"); }
    const char*     Title()             { return GetTextValue("Title"); }
    const char*     Subject()           { return GetTextValue("Subject"); }
    const char*     Keywords()          { return GetTextValue("Keywords"); }
    void            SetCreationDate(pdf_date value);
    void            SetModDate(pdf_date value);
    void            SetAuthor(const char* value, 
                        PdfEncodingDef* encoding = NULL);
    void            SetCreator(const char* value, 
                        PdfEncodingDef* encoding = NULL);
    void            SetProducer(const char* value, 
                        PdfEncodingDef* encoding = NULL);
    void            SetTitle(const char* value, 
                        PdfEncodingDef* encoding = NULL);
    void            SetSubject(const char* value, 
                        PdfEncodingDef* encoding = NULL);
    void            SetKeywords(const char* value, 
                        PdfEncodingDef* encoding = NULL);
    void            SetAuthorMb(const char* value, PdfCMap* cmap)
                                { SetMbText("Author", value, cmap); }
    void            SetCreatorMb(const char* value, PdfCMap* cmap)
                                { SetMbText("Creator", value, cmap); }
    void            SetProducerMb(const char* value, PdfCMap* cmap)
                                { SetMbText("Producer", value, cmap); }
    void            SetTitleMb(const char* value, PdfCMap* cmap)
                                { SetMbText("Title", value, cmap); }
    void            SetSubjectMb(const char* value, PdfCMap* cmap)
                                { SetMbText("Subject", value, cmap); }
    void            SetKeywordsMb(const char* value, PdfCMap* cmap)
                                { SetMbText("Keywords", value, cmap); }
protected:
private:
    pdf_date        StrToPdfDate(const char* value);
    bool            PdfDateToStr(const pdf_date date, char* value, int length);
    void            SetTextAsUnicode(const char* key, const char* value,
                        PdfEncodingDef* encoding);
    void            SetMbText(const char* key, const char* value, 
                        PdfCMap* cmap);
};

/*----------------------------------------------------------------------------*/
/*----- PdfEncryptDict class -------------------------------------------------*/

#ifdef USE_ENCRYPTION
class PdfEncryptDict : public PdfDictionary
{
public:
                    PdfEncryptDict(PdfXref* xref, PdfInfo* info = NULL);
    void            Init();
    void            SetPassword(const char* owner_passwd, 
                        const char* user_passwd);
    void            SetPermission(int value)
                        { fPermission = value | PDF_PERMISSION_PAD; }
    void            EncryptPassword();
    int             Permission()    { return fPermission; }
    const unsigned char*
                    ID()    { return fID; }
    const unsigned char*
                    FileKey()       { return fFileKey; }
private:
    void            CreateID();
    void            CreateUserKey();
    void            CreateOwnerKey();
    void            CreateFileKey();
    void            PadOrTruncatePasswd(const char* passwd, unsigned char* out);
    unsigned char   fOwnerPasswdValue[PDF_PASSWD_LEN];
    unsigned char   fUserPasswdValue[PDF_PASSWD_LEN];
    unsigned char   fOwnerPasswd[PDF_PASSWD_LEN];
    unsigned char   fUserPasswd[PDF_PASSWD_LEN];
    unsigned char   fFileKey[PDF_MD5_KEY_LEN];
    unsigned char   fID[PDF_ID_LEN];
    PdfInfo*        fInfo;
    int             fPermission;
};
#else
typedef void PdfEncryptDict;
#endif

/*----------------------------------------------------------------------------*/
/*----- PdfOutlineBase class -------------------------------------------------*/

class PdfOutlineItem;

class PdfOutlineBase : public PdfDictionary
{
    friend class PdfOutlineItem;
public:
                        PdfOutlineBase(PdfXref *xref);
        PdfOutlineItem* First()
                        { return (PdfOutlineItem*)GetValue("First"); }
        PdfOutlineItem* Last()
                        { return (PdfOutlineItem*)GetValue("Last"); }
        bool            HasChild()    { return ChildCount() > 0; }
virtual int             ChildCount();
        bool            Opened()      { return fOpened; }
        void            SetOpened(bool value)   { fOpened = value; }
protected:
        void            AddChild(PdfOutlineItem* item);
        void            InternalWriteStream(PdfStreamBase* out,
                            PdfEncryptor* e);
        bool            fOpened;
};

typedef PdfOutlineBase  PdfOutline;

/*----------------------------------------------------------------------------*/
/*----- PdfOutlineRoot class -------------------------------------------------*/

class PdfOutlineRoot : public PdfOutline
{
public:
                        PdfOutlineRoot(PdfXref *xref);
};

/*----------------------------------------------------------------------------*/
/*----- PdfOutlineItem class -------------------------------------------------*/

class PdfOutlineItem : public PdfOutline
{
public:
                        PdfOutlineItem(PdfOutline *parent);
        const char*     Title();
        void            SetTitle(const char* title, 
                            PdfEncodingDef* encoding = NULL);
        void            SetTitleMb(const char* title, PdfCMap* cmap);
        void            SetDestination(PdfDestination* dst);
        PdfOutlineItem* Prev()
                        { return (PdfOutlineItem*)GetValue("Prev"); }
        PdfOutlineItem* Next()
                        { return (PdfOutlineItem*)GetValue("Next"); }
        PdfOutline*     Parent()
                        { return (PdfOutline*)GetValue("Parent"); }
private:
        PdfText*        fTitle;
};

/*----------------------------------------------------------------------------*/
/*----- PdfXObjectBase class -------------------------------------------------*/

class PdfXObjectBase : public PdfStream
{
    friend class PdfXObjectMgr;
public:
                        PdfXObjectBase(PdfDoc* doc);
                        ~PdfXObjectBase();
        const char*     Name()    { return fName; }
virtual bool            IsValidObject() = 0;
protected:
        void            SetName(const char* name);
private:
        char*           fName;
        bool            fValid;
};

/*----------------------------------------------------------------------------*/
/*----- PdfXObjectMgr class --------------------------------------------------*/

class PdfXObjectMgr
{
public:
                        PdfXObjectMgr(PdfXref* xref);
                        ~PdfXObjectMgr();
        void            RegisterXObject(PdfXObject* xobject, const char* name);
        PdfXObject*     GetXObject(const char* name);
        PdfXObject*     GetXObject(unsigned int index)
                            { return (PdfXObject*)fList->ItemAt(index); }
        unsigned int    CountXObjects()  { return fList->CountItems(); }
private:
        PdfList*        fList;
        PdfXref*        fXref;
};

/*----------------------------------------------------------------------------*/
/*----- PdfImage class -------------------------------------------------------*/

class PdfImage : public PdfXObject
{
public:
                        PdfImage(PdfDoc* doc);
        double           Width()              { return fWidth; }
        double           Height()             { return fHeight; }
        pdf_color_space ColorSpace()         { return fColorSpace; }
        unsigned int    BitsPerComponent()   { return fBitsPerComponent; }
protected:
        void            InternalWriteStream(PdfStreamBase* out,
                            PdfEncryptor* e);
        unsigned int    fWidth;
        unsigned int    fHeight;
        unsigned int    fBitsPerComponent;
        pdf_color_space fColorSpace;
        pdf_pal_color*  fPallet;
        unsigned int    fNumPallet;
};

/*----------------------------------------------------------------------------*/
/*------ utilities -----------------------------------------------------------*/

const char* PdfGetColorSpaceName(pdf_color_space cs);

/*----------------------------------------------------------------------------*/
/*----- PdfHeader class ------------------------------------------------------*/

class PdfHeader
{
    friend class    PdfDoc;
protected:
    void            WriteToStream(PdfStreamBase* out);
private:
};

/*----------------------------------------------------------------------------*/
/*----- PdfTrailer class -----------------------------------------------------*/

class PdfTrailer
{
    friend class    PdfDoc;
public:
                    PdfTrailer(PdfXref *xref);
                    ~PdfTrailer();
    void            SetEncryptDict(PdfEncryptDict* encrypt);
protected:
    void            WriteToStream(PdfStreamBase* out);
    void            SetXrefAddr(unsigned int addr);
    void            SetXrefSize(unsigned int size);
    void            SetRoot(PdfCatalog* root);
    void            SetInfo(PdfInfo* info);
private:
    void            SetID(const unsigned char* id);
    PdfDictionary*  fAttributes;
    PdfEncryptDict* fEncryptDict;
    unsigned int    fXrefAddr;
    PdfXref*        fXref;
    unsigned char   fID1[PDF_ID_LEN];
    unsigned char   fID2[PDF_ID_LEN];
};

/*----------------------------------------------------------------------------*/
/*----- PdfDoc class ---------------------------------------------------------*/

class PdfDoc
{
    friend class PdfPages;
public:
                        PdfDoc();
                        ~PdfDoc();
        void            NewDoc();
        void            FreeDoc(bool free_all_objects = true);
        void            WriteToStream(PdfStreamBase* out);
        void            WriteToFile(const char* filename);
        void            AddType1Font(PdfType1FontDef* fontdef, 
                            const char* name = NULL, 
                            PdfEncodingDef* encoding = NULL);

        void            AddType0Font(PdfCIDFontDef* fontdef, 
                            const char* name = NULL, 
                            PdfCMap* cmap = NULL);

        void            AddXObject(PdfXObject* xobject, 
                            const char* name = NULL);
        void            RegisterObject(PdfAutoPtrObject* obj);
        PdfPage*        AddPage();
        PdfCatalog*     Catalog()           { return fCatalog; }
        PdfInfo*        Info();
        void            SetPassword(const char* owner_passwd,
                            const char* user_passwd);
        void            SetPermission(int permission);
        PdfXref*        Xref()              { return fXref; }
        PdfOutlineRoot* Outlines()          { return Catalog()->Outlines(); }
        PdfPages*       RootPages()         { return fRootPages; }
        PdfPages*       CurrentPages()      { return fCurrentPages; }
        PdfPage*        CurrentPage()       { return fCurrentPage; }
        PdfXObjectMgr*  XObjectMgr()        { return fXObjectMgr; }
        PdfFontMgr*     FontMgr()           { return fFontMgr; }
        bool            HasDoc()            { return fHasDoc; }
        int             LastError()         { return fError; }
        void            SetError(int err);
private:
        void            Init();
        void            SetCurrentPages(PdfPages* pages)
                                            { fCurrentPages = pages; }
        void            SetCurrentPage(PdfPage* page)
                                            { fCurrentPage = page; }
        PdfEncryptDict* GetEncryptDict()
                            { return fTrailer->fEncryptDict; }
        bool            IsEncrypted()   
                            { return (GetEncryptDict() != NULL); }
        PdfXref*        fXref;
        PdfHeader*      fHeader;
        PdfTrailer*     fTrailer;
        PdfCatalog*     fCatalog;
        PdfInfo*        fInfo;
        PdfFontMgr*     fFontMgr;
        PdfXObjectMgr*  fXObjectMgr;
        PdfAutoPtrMgr*  fAutoPtrMgr;
        PdfPages*       fRootPages;
        PdfPages*       fCurrentPages;
        PdfPage*        fCurrentPage;
        bool            fHasDoc;
        int             fError;
};

/*----------------------------------------------------------------------------*/
/*----- PdfPngImage class ----------------------------------------------------*/

#ifndef NOPNG

class PdfPngImage : public PdfImage
{
public:
                        PdfPngImage(PdfDoc* doc);
                        ~PdfPngImage();
        void            LoadFromFile(const char* filename);
        void            FreeImage();
        bool            IsValidObject()    { return fHasImage; }
private:
        FILE*           OpenImageFile(const char* filename);
        void            CreatePallet(png_structp png_ptr, png_infop info_ptr);
        bool            fHasImage;
};

#endif /* NOPNG */

/*----------------------------------------------------------------------------*/
/*----- PdfJpegImage class ---------------------------------------------------*/

#ifndef NOJPEG

class PdfJpegImage : public PdfImage
{
public:
                        PdfJpegImage(PdfDoc* doc);
                        ~PdfJpegImage();
        void            LoadFromFile(const char* filename);
        void            FreeImage();
        bool            IsValidObject()    { return fHasImage; }
private:
        bool            fHasImage;
};

void
PdfJpegErrorExit(j_common_ptr cinfo);

#endif /* NOJPEG */


/*---------------------------------------------------------------------------*/
/*----- Utility routines ----------------------------------------------------*/

extern "C" {
#endif /* __cplusplus */

pdf_rgb_color PdfRGBColor(double r, double g, double b);

pdf_point PdfPoint(double x, double y);

pdf_rect PdfRect(double left, double bottom, double right, double top);

pdf_box PdfBox(int left, int bottom, int right, int top);

pdf_text_matrix PdfTextMatrix(double a, double b, double c, double d,
        double x, double y);

#ifdef __cplusplus
}
#endif

#endif /* _LIB_HARU_H */

