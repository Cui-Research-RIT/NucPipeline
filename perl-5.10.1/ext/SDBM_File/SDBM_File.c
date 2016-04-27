/*
 * This file was generated automatically by ExtUtils::ParseXS version 2.2002 from the
 * contents of SDBM_File.xs. Do not edit this file, edit SDBM_File.xs instead.
 *
 *	ANY CHANGES MADE HERE WILL BE LOST! 
 *
 */

#line 1 "SDBM_File.xs"
#define PERL_NO_GET_CONTEXT
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "sdbm/sdbm.h"

typedef struct {
	DBM * 	dbp ;
	SV *    filter_fetch_key ;
	SV *    filter_store_key ;
	SV *    filter_fetch_value ;
	SV *    filter_store_value ;
	int     filtering ;
	} SDBM_File_type;

typedef SDBM_File_type * SDBM_File ;
typedef datum datum_key ;
typedef datum datum_value ;

#define sdbm_TIEHASH(dbtype,filename,flags,mode) sdbm_open(filename,flags,mode)
#define sdbm_FETCH(db,key)			sdbm_fetch(db->dbp,key)
#define sdbm_STORE(db,key,value,flags)		sdbm_store(db->dbp,key,value,flags)
#define sdbm_DELETE(db,key)			sdbm_delete(db->dbp,key)
#define sdbm_EXISTS(db,key)			sdbm_exists(db->dbp,key)
#define sdbm_FIRSTKEY(db)			sdbm_firstkey(db->dbp)
#define sdbm_NEXTKEY(db,key)			sdbm_nextkey(db->dbp)


#line 39 "SDBM_File.c"
#ifndef PERL_UNUSED_VAR
#  define PERL_UNUSED_VAR(var) if (0) var = var
#endif

#ifndef PERL_ARGS_ASSERT_CROAK_XS_USAGE
#define PERL_ARGS_ASSERT_CROAK_XS_USAGE assert(cv); assert(params)

/* prototype to pass -Wmissing-prototypes */
STATIC void
S_croak_xs_usage(pTHX_ const CV *const cv, const char *const params);

STATIC void
S_croak_xs_usage(pTHX_ const CV *const cv, const char *const params)
{
    const GV *const gv = CvGV(cv);

    PERL_ARGS_ASSERT_CROAK_XS_USAGE;

    if (gv) {
        const char *const gvname = GvNAME(gv);
        const HV *const stash = GvSTASH(gv);
        const char *const hvname = stash ? HvNAME(stash) : NULL;

        if (hvname)
            Perl_croak(aTHX_ "Usage: %s::%s(%s)", hvname, gvname, params);
        else
            Perl_croak(aTHX_ "Usage: %s(%s)", gvname, params);
    } else {
        /* Pants. I don't think that it should be possible to get here. */
        Perl_croak(aTHX_ "Usage: CODE(0x%"UVxf")(%s)", PTR2UV(cv), params);
    }
}
#undef  PERL_ARGS_ASSERT_CROAK_XS_USAGE

#ifdef PERL_IMPLICIT_CONTEXT
#define croak_xs_usage(a,b)	S_croak_xs_usage(aTHX_ a,b)
#else
#define croak_xs_usage		S_croak_xs_usage
#endif

#endif

#line 82 "SDBM_File.c"

XS(XS_SDBM_File_TIEHASH); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_TIEHASH)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 4)
       croak_xs_usage(cv,  "dbtype, filename, flags, mode");
    {
	char *	dbtype = (char *)SvPV_nolen(ST(0));
	char *	filename = (char *)SvPV_nolen(ST(1));
	int	flags = (int)SvIV(ST(2));
	int	mode = (int)SvIV(ST(3));
	SDBM_File	RETVAL;
#line 38 "SDBM_File.xs"
	{
	    DBM * 	dbp ;

	    RETVAL = NULL ;
	    if ((dbp = sdbm_open(filename,flags,mode))) {
	        RETVAL = (SDBM_File)safemalloc(sizeof(SDBM_File_type)) ;
    	        Zero(RETVAL, 1, SDBM_File_type) ;
		RETVAL->dbp = dbp ;
	    }

	}
#line 112 "SDBM_File.c"
	ST(0) = sv_newmortal();
        sv_setref_pv(ST(0), dbtype, (void*)RETVAL);
    }
    XSRETURN(1);
}


XS(XS_SDBM_File_DESTROY); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_DESTROY)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 1)
       croak_xs_usage(cv,  "db");
    {
	SDBM_File	db;

	if (SvROK(ST(0))) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not a reference",
			"SDBM_File::DESTROY",
			"db");
#line 56 "SDBM_File.xs"
	if (db) {
	    sdbm_close(db->dbp);
	    if (db->filter_fetch_key)
		SvREFCNT_dec(db->filter_fetch_key) ;
	    if (db->filter_store_key)
		SvREFCNT_dec(db->filter_store_key) ;
	    if (db->filter_fetch_value)
		SvREFCNT_dec(db->filter_fetch_value) ;
	    if (db->filter_store_value)
		SvREFCNT_dec(db->filter_store_value) ;
	    safefree(db) ;
	}
#line 154 "SDBM_File.c"
    }
    XSRETURN_EMPTY;
}


XS(XS_SDBM_File_FETCH); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_FETCH)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 2)
       croak_xs_usage(cv,  "db, key");
    {
	SDBM_File	db;
	datum_key	key;
	datum_value	RETVAL;

	if (sv_derived_from(ST(0), "SDBM_File")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not of type %s",
			"SDBM_File::FETCH",
			"db", "SDBM_File");

	{
	    STRLEN len;
	    DBM_ckFilter(ST(1), filter_store_key, "filter_store_key");
	    key.dptr = SvPVbyte(ST(1), len);
	    key.dsize = (int)len;
	};

	RETVAL = sdbm_FETCH(db, key);
	ST(0) = sv_newmortal();
	sv_setpvn(ST(0), RETVAL.dptr, RETVAL.dsize);
	DBM_ckFilter(ST(0), filter_fetch_value,"filter_fetch_value");
    }
    XSRETURN(1);
}


XS(XS_SDBM_File_STORE); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_STORE)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items < 3 || items > 4)
       croak_xs_usage(cv,  "db, key, value, flags = DBM_REPLACE");
    {
	SDBM_File	db;
	datum_key	key;
	datum_value	value;
	int	flags;
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(0), "SDBM_File")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not of type %s",
			"SDBM_File::STORE",
			"db", "SDBM_File");

	{
	    STRLEN len;
	    DBM_ckFilter(ST(1), filter_store_key, "filter_store_key");
	    key.dptr = SvPVbyte(ST(1), len);
	    key.dsize = (int)len;
	};

        DBM_ckFilter(ST(2), filter_store_value, "filter_store_value");
	if (SvOK(ST(2))) {
	    STRLEN len;
	    value.dptr = SvPVbyte(ST(2), len);
	    value.dsize = (int)len;
	}
	else {
	    value.dptr = "";
	    value.dsize = 0;
	};

	if (items < 4)
	    flags = DBM_REPLACE;
	else {
	    flags = (int)SvIV(ST(3));
	}

	RETVAL = sdbm_STORE(db, key, value, flags);
	XSprePUSH; PUSHi((IV)RETVAL);
#line 81 "SDBM_File.xs"
	if (RETVAL) {
	    if (RETVAL < 0 && errno == EPERM)
		croak("No write permission to sdbm file");
	    croak("sdbm store returned %d, errno %d, key \"%s\"",
			RETVAL,errno,key.dptr);
	    sdbm_clearerr(db->dbp);
	}
#line 261 "SDBM_File.c"
    }
    XSRETURN(1);
}


XS(XS_SDBM_File_DELETE); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_DELETE)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 2)
       croak_xs_usage(cv,  "db, key");
    {
	SDBM_File	db;
	datum_key	key;
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(0), "SDBM_File")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not of type %s",
			"SDBM_File::DELETE",
			"db", "SDBM_File");

	{
	    STRLEN len;
	    DBM_ckFilter(ST(1), filter_store_key, "filter_store_key");
	    key.dptr = SvPVbyte(ST(1), len);
	    key.dsize = (int)len;
	};

	RETVAL = sdbm_DELETE(db, key);
	XSprePUSH; PUSHi((IV)RETVAL);
    }
    XSRETURN(1);
}


XS(XS_SDBM_File_EXISTS); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_EXISTS)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 2)
       croak_xs_usage(cv,  "db, key");
    {
	SDBM_File	db;
	datum_key	key;
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(0), "SDBM_File")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not of type %s",
			"SDBM_File::EXISTS",
			"db", "SDBM_File");

	{
	    STRLEN len;
	    DBM_ckFilter(ST(1), filter_store_key, "filter_store_key");
	    key.dptr = SvPVbyte(ST(1), len);
	    key.dsize = (int)len;
	};

	RETVAL = sdbm_EXISTS(db, key);
	XSprePUSH; PUSHi((IV)RETVAL);
    }
    XSRETURN(1);
}


XS(XS_SDBM_File_FIRSTKEY); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_FIRSTKEY)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 1)
       croak_xs_usage(cv,  "db");
    {
	SDBM_File	db;
	datum_key	RETVAL;

	if (sv_derived_from(ST(0), "SDBM_File")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not of type %s",
			"SDBM_File::FIRSTKEY",
			"db", "SDBM_File");

	RETVAL = sdbm_FIRSTKEY(db);
	ST(0) = sv_newmortal();
	sv_setpvn(ST(0), RETVAL.dptr, RETVAL.dsize);
	DBM_ckFilter(ST(0), filter_fetch_key,"filter_fetch_key");
    }
    XSRETURN(1);
}


XS(XS_SDBM_File_NEXTKEY); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_NEXTKEY)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 2)
       croak_xs_usage(cv,  "db, key");
    {
	SDBM_File	db;
	datum_key	key;
	datum_key	RETVAL;

	if (sv_derived_from(ST(0), "SDBM_File")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not of type %s",
			"SDBM_File::NEXTKEY",
			"db", "SDBM_File");

	{
	    STRLEN len;
	    DBM_ckFilter(ST(1), filter_store_key, "filter_store_key");
	    key.dptr = SvPVbyte(ST(1), len);
	    key.dsize = (int)len;
	};

	RETVAL = sdbm_NEXTKEY(db, key);
	ST(0) = sv_newmortal();
	sv_setpvn(ST(0), RETVAL.dptr, RETVAL.dsize);
	DBM_ckFilter(ST(0), filter_fetch_key,"filter_fetch_key");
    }
    XSRETURN(1);
}


XS(XS_SDBM_File_error); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_error)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 1)
       croak_xs_usage(cv,  "db");
    {
	SDBM_File	db;
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(0), "SDBM_File")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not of type %s",
			"SDBM_File::error",
			"db", "SDBM_File");
#line 112 "SDBM_File.xs"
	RETVAL = sdbm_error(db->dbp) ;
#line 442 "SDBM_File.c"
	XSprePUSH; PUSHi((IV)RETVAL);
    }
    XSRETURN(1);
}


XS(XS_SDBM_File_clearerr); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_clearerr)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 1)
       croak_xs_usage(cv,  "db");
    {
	SDBM_File	db;
	int	RETVAL;
	dXSTARG;

	if (sv_derived_from(ST(0), "SDBM_File")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not of type %s",
			"SDBM_File::clearerr",
			"db", "SDBM_File");
#line 120 "SDBM_File.xs"
	RETVAL = sdbm_clearerr(db->dbp) ;
#line 474 "SDBM_File.c"
	XSprePUSH; PUSHi((IV)RETVAL);
    }
    XSRETURN(1);
}


XS(XS_SDBM_File_filter_fetch_key); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_filter_fetch_key)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 2)
       croak_xs_usage(cv,  "db, code");
    {
	SDBM_File	db;
	SV *	code = ST(1);
	SV *	RETVAL = &PL_sv_undef ;

	if (sv_derived_from(ST(0), "SDBM_File")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not of type %s",
			"SDBM_File::filter_fetch_key",
			"db", "SDBM_File");
#line 131 "SDBM_File.xs"
	    DBM_setFilter(db->filter_fetch_key, code) ;
#line 506 "SDBM_File.c"
    }
    XSRETURN(1);
}


XS(XS_SDBM_File_filter_store_key); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_filter_store_key)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 2)
       croak_xs_usage(cv,  "db, code");
    {
	SDBM_File	db;
	SV *	code = ST(1);
	SV *	RETVAL =  &PL_sv_undef ;

	if (sv_derived_from(ST(0), "SDBM_File")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not of type %s",
			"SDBM_File::filter_store_key",
			"db", "SDBM_File");
#line 139 "SDBM_File.xs"
	    DBM_setFilter(db->filter_store_key, code) ;
#line 537 "SDBM_File.c"
    }
    XSRETURN(1);
}


XS(XS_SDBM_File_filter_fetch_value); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_filter_fetch_value)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 2)
       croak_xs_usage(cv,  "db, code");
    {
	SDBM_File	db;
	SV *	code = ST(1);
	SV *	RETVAL =  &PL_sv_undef ;

	if (sv_derived_from(ST(0), "SDBM_File")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not of type %s",
			"SDBM_File::filter_fetch_value",
			"db", "SDBM_File");
#line 147 "SDBM_File.xs"
	    DBM_setFilter(db->filter_fetch_value, code) ;
#line 568 "SDBM_File.c"
    }
    XSRETURN(1);
}


XS(XS_SDBM_File_filter_store_value); /* prototype to pass -Wmissing-prototypes */
XS(XS_SDBM_File_filter_store_value)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    if (items != 2)
       croak_xs_usage(cv,  "db, code");
    {
	SDBM_File	db;
	SV *	code = ST(1);
	SV *	RETVAL =  &PL_sv_undef ;

	if (sv_derived_from(ST(0), "SDBM_File")) {
	    IV tmp = SvIV((SV*)SvRV(ST(0)));
	    db = INT2PTR(SDBM_File,tmp);
	}
	else
	    Perl_croak(aTHX_ "%s: %s is not of type %s",
			"SDBM_File::filter_store_value",
			"db", "SDBM_File");
#line 155 "SDBM_File.xs"
	    DBM_setFilter(db->filter_store_value, code) ;
#line 599 "SDBM_File.c"
    }
    XSRETURN(1);
}

#ifdef __cplusplus
extern "C"
#endif
XS(boot_SDBM_File); /* prototype to pass -Wmissing-prototypes */
XS(boot_SDBM_File)
{
#ifdef dVAR
    dVAR; dXSARGS;
#else
    dXSARGS;
#endif
    const char* file = __FILE__;

    PERL_UNUSED_VAR(cv); /* -W */
    PERL_UNUSED_VAR(items); /* -W */
    XS_VERSION_BOOTCHECK ;

        newXS("SDBM_File::TIEHASH", XS_SDBM_File_TIEHASH, file);
        newXS("SDBM_File::DESTROY", XS_SDBM_File_DESTROY, file);
        newXS("SDBM_File::FETCH", XS_SDBM_File_FETCH, file);
        newXS("SDBM_File::STORE", XS_SDBM_File_STORE, file);
        newXS("SDBM_File::DELETE", XS_SDBM_File_DELETE, file);
        newXS("SDBM_File::EXISTS", XS_SDBM_File_EXISTS, file);
        newXS("SDBM_File::FIRSTKEY", XS_SDBM_File_FIRSTKEY, file);
        newXS("SDBM_File::NEXTKEY", XS_SDBM_File_NEXTKEY, file);
        newXS("SDBM_File::error", XS_SDBM_File_error, file);
        newXS("SDBM_File::clearerr", XS_SDBM_File_clearerr, file);
        newXS("SDBM_File::filter_fetch_key", XS_SDBM_File_filter_fetch_key, file);
        newXS("SDBM_File::filter_store_key", XS_SDBM_File_filter_store_key, file);
        newXS("SDBM_File::filter_fetch_value", XS_SDBM_File_filter_fetch_value, file);
        newXS("SDBM_File::filter_store_value", XS_SDBM_File_filter_store_value, file);
    if (PL_unitcheckav)
         call_list(PL_scopestack_ix, PL_unitcheckav);
    XSRETURN_YES;
}

