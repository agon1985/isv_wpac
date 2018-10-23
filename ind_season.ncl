;******************************************************************************
; File name: ind_season.ncl
; Author:    Alex Gonzalez
; Date:      December 2017
; Contact:   alex.o.gonzalez@jpl.nasa.gov
;******************************************************************************

undef("ind_season")

function ind_season( inData[*]:numeric, monthStart:integer, \\
                    monthEnd:integer, yearStart:integer, yearEnd:integer )
;******************************************************************************
; Apply time filtering of data to a season of interest by computing the indices
; associated with the season. The season can be as short as one month and as 
; long as 12 months. Also, the season can span two years, e.g., December 
; through February, and this is represented by cutting off two partial seasons
; associated with the first and last year. 
; 
;
; Input Variables:
;   inData: the data must be temporal data and can be of any numeric type
;   monthStart: integer representing first month of season
;   monthEnd:   integer representing last month of season
;   yearStart:  integer representing first year of season
;   yearEnd:    integer representing last year of season
;
; Return Value:
;   retVal[*][*]: the indices associated with the season
;******************************************************************************

local mndiff, all_months, time, t3d, lmon, lyr, id_seas, retVal

begin
; e.g., if season is Feb–Jan or Mar–Feb, it's considered whole year
; and mndiff = -1, or if season is Jan–Dec then mndiff = 11
mndiff = (/-1,11/)
if any(monthEnd-monthStart .eq. mndiff)
    all_months = True
    print("FYI, the season you chose is the entire year...")
else
    all_months = False
end if

time = inData
printVarSummary(time)
; make a long vector with years and months for each time 
t3d = cd_calendar(time,0)
lyr  = t3d(:,0)
lmon = t3d(:,1)
delete([/t3d,time/])

; find indices associated with chosen years and months
if (all_months .eq. False)
   if (monthEnd .gt. monthStart)      ; season spans one year
      id_seas = ind( (lyr.ge.yearStart.and.lyr.le.yearEnd) .and. (lmon.ge.monthStart.and.lmon.le.monthEnd) )
   else if (monthEnd .lt. monthStart) ; season spans two years
      ; the first and last years are NOT considered because they both 
      ; only have part of a season
      id_seas = ind( (lyr.ge.yearStart.and.lyr.lt.yearEnd.and.lmon.ge.monthStart) \
                     .or. (lyr.gt.yearStart.and.lyr.le.yearEnd.and.lmon.le.monthEnd) )
      ; split into two conditional statements, for first and second part of season
   else ; monthEnd .eq. monthStart ; only one month
      id_seas = ind( (lyr.ge.yearStart.and.lyr.le.yearEnd) .and. (lmon.eq.monthStart) )
   end if
   end if
else ;(all_months .eq. True) so just choose years of interest
   id_seas = ind(lyr.ge.yearStart.and.lyr.le.yearEnd)
end if
delete([/lyr,lmon/])
;printVarSummary(id_seas)
retVal = id_seas

  ; return the filtered data
  return( retVal )
end
